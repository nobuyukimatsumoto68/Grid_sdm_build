/*
 * wflow_q_h5_claude.cc
 *
 * Gluonic topological charge from the Wilson (gradient) flow, one h5 file per gauge configuration.
 * Flow logic as Grid's HMC/ComputeWilsonFlow.cc, but with the library-default flow action
 * WilsonGaugeAction(Nc) (correct normalization for SU(4); ComputeWilsonFlow's flow_type "Wilson" sets
 * WilsonGaugeAction(3.0), i.e. flow time x 3/4 for Nc = 4) and h5 output instead of stdout only.
 * Plan: Grid/scripts/wflow_topo_q_impl_plan_claude.md (round 2).
 *
 * Methods:
 *   Wilson flow, RK3 integrator: M. Luescher, JHEP 08 (2010) 071, arXiv:1006.4518 (Grid WilsonFlow::evolve_step).
 *   Q clover: WilsonLoops::TopologicalCharge; Q 5Li: P. de Forcrand, M. Garcia Perez, I.-O. Stamatescu,
 *   NPB 499 (1997) 409, hep-lat/9701012 (WilsonLoops::TopologicalCharge5Li).
 *
 * Per config <outdir>/wflow.<conf>.h5 (Hdf5Writer), arrays over the measurement steps (step 0 and every MEAS):
 *   t, t2E_plaq (= t^2 <E> plaquette form), t2E_clover, Q_clover, Q5Li, plus plaq0 (t = 0 plaquette).
 * Written to <file>.tmp and renamed only when complete (a hang never leaves a partial wflow.<conf>.h5);
 * configs with an existing wflow.<conf>.h5 are skipped (rerun = resume); missing gauge files are skipped.
 * Config selection: conf = start, start+step, ..., <= end, and only every nparts-th of those (index % nparts ==
 * part), so several single-rank processes on one node share a range.
 *
 * Usage (single rank per process):
 *   wflow_q_h5_claude --dir <latdir> --prefix <conf_..._lat> --start S --end E --step K --outdir <dir>
 *                     [--part p --nparts n] --grid 24.24.24.48 --mpi 1.1.1.1 ...
 * Env: WFLOW_STEPS (600), WFLOW_STEP_SIZE (0.01), WFLOW_MEAS (20), WFLOW_DEADLINE_EPOCH (graceful stop, 0 = off),
 *      WFLOW_TPT_SECONDS (bootstrap per-config wall estimate).
 */
#include <Grid/Grid.h>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <string>
#include <vector>

using namespace Grid;

static std::string arg_or(int argc, char** argv, const std::string& opt, const std::string& def)
{
  if (GridCmdOptionExists(argv, argv + argc, opt)) return GridCmdOptionPayload(argv, argv + argc, opt);
  return def;
}

static double env_or(const char* name, double def)
{
  const char* e = std::getenv(name);
  if (e == nullptr) return def;
  return std::atof(e);
}

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  const std::string dir = arg_or(argc, argv, "--dir", "");
  const std::string prefix = arg_or(argc, argv, "--prefix", "");
  const std::string outdir = arg_or(argc, argv, "--outdir", "");
  const int start = std::atoi(arg_or(argc, argv, "--start", "0").c_str());
  const int end = std::atoi(arg_or(argc, argv, "--end", "-1").c_str());
  const int cstep = std::atoi(arg_or(argc, argv, "--step", "1").c_str());
  const int part = std::atoi(arg_or(argc, argv, "--part", "0").c_str());
  const int nparts = std::atoi(arg_or(argc, argv, "--nparts", "1").c_str());
  GRID_ASSERT(!dir.empty() && !prefix.empty() && !outdir.empty() && cstep > 0 && nparts > 0);

  const int nsteps = (int)env_or("WFLOW_STEPS", 600);
  const double eps = env_or("WFLOW_STEP_SIZE", 0.01);
  const int meas = (int)env_or("WFLOW_MEAS", 20);
  const long deadline = (long)env_or("WFLOW_DEADLINE_EPOCH", 0);
  double max_dur = env_or("WFLOW_TPT_SECONDS", 0.0);
  const double margin = 1.2;

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                        GridDefaultSimd(Nd, vComplex::Nsimd()),
                                                        GridDefaultMpi());
  LatticeGaugeField Umu(UGrid);
  LatticeGaugeField Uflow(UGrid);
  std::filesystem::create_directories(outdir);

  std::cout << GridLogMessage << "# wflow_q_h5: dir=" << dir << " prefix=" << prefix << " range " << start << ".."
            << end << " step " << cstep << " part " << part << "/" << nparts << " outdir=" << outdir
            << " flow steps=" << nsteps << " eps=" << eps << " meas=" << meas << std::endl;

  int index = 0;
  for (int conf = start; conf <= end; conf += cstep) {
    const int my = index % nparts;
    index++;
    if (my != part) continue;

    const std::string cfg = dir + "/" + prefix + "." + std::to_string(conf);
    const std::string out = outdir + "/wflow." + std::to_string(conf) + ".h5";
    if (std::filesystem::exists(out)) {
      std::cout << GridLogMessage << "# skipping conf " << conf << " (output exists)" << std::endl;
      continue;
    }
    if (!std::filesystem::exists(cfg)) {
      std::cout << GridLogMessage << "# skipping conf " << conf << " (no gauge file " << cfg << ")" << std::endl;
      continue;
    }
    if (deadline > 0 && max_dur > 0.0) {
      uint64_t stop = 0;
      if (UGrid->IsBoss()) {
        if ((double)std::time(nullptr) + margin * max_dur > (double)deadline) stop = 1;
      }
      UGrid->GlobalSum(stop);
      if (stop) {
        std::cout << GridLogMessage << "# blocker: stopping before conf " << conf << std::endl;
        break;
      }
    }
    const long t0 = (long)std::time(nullptr);

    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, cfg);
    const RealD plaq0 = WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);

    std::vector<RealD> tv;
    std::vector<RealD> t2Ep;
    std::vector<RealD> t2Ec;
    std::vector<RealD> qc;
    std::vector<RealD> q5;
    // library-default flow action WilsonGaugeAction(Nc) (no setGaugeAction call)
    WilsonFlow<PeriodicGimplR> WF(eps, nsteps, meas);
    WF.resetActions();
    // the measurement callback must capture the output vectors (std::function interface of WilsonFlow)
    WF.addMeasurement(meas, [&tv, &t2Ep, &t2Ec, &qc, &q5](int step, RealD t, const PeriodicGimplR::GaugeField& U) {
      tv.push_back(t);
      t2Ep.push_back(WilsonFlow<PeriodicGimplR>::energyDensityPlaquette(t, U));
      t2Ec.push_back(WilsonFlow<PeriodicGimplR>::energyDensityCloverleaf(t, U));
      qc.push_back(WilsonLoops<PeriodicGimplR>::TopologicalCharge(U));
      q5.push_back(WilsonLoops<PeriodicGimplR>::TopologicalCharge5Li(U));
    });
    WF.smear(Uflow, Umu);

    if (UGrid->IsBoss()) {
      const std::string tmp = out + ".tmp";
      {
        Hdf5Writer WR(tmp);
        write(WR, "t", tv);
        write(WR, "t2E_plaq", t2Ep);
        write(WR, "t2E_clover", t2Ec);
        write(WR, "Q_clover", qc);
        write(WR, "Q5Li", q5);
        write(WR, "plaq0", plaq0);
        write(WR, "conf", conf);
      }
      std::filesystem::rename(tmp, out);
    }

    uint64_t dur = 0;
    if (UGrid->IsBoss()) dur = (uint64_t)((long)std::time(nullptr) - t0);
    UGrid->GlobalSum(dur);
    if ((double)dur > max_dur) max_dur = (double)dur;
    std::cout << GridLogMessage << "# conf " << conf << " took " << dur << "s  Q5Li(t=" << tv.back() << ")="
              << q5.back() << "  Q_clover=" << qc.back() << std::endl;
  }
  std::cout << GridLogMessage << "# wflow_q_h5 done." << std::endl;
  Grid_finalize();
  return 0;
}
