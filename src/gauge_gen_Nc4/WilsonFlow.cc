#include <Grid/Grid.h>
// 082224 sungwoo:
// copied from /usr/WS2/lsd/matsumoto5/su3/Grid/examples/
// originally dated Jun 20 2024

using namespace std;
using namespace Grid;


class FlowResult: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(FlowResult,
				  std::string, description,
				  std::vector<RealD>, data);
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                        GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                        GridDefaultMpi());

  LatticeGaugeField U(UGrid);
  std::string config;
  int interval = 2;
  int nsteps = 200;
  double step_size = 0.01;
  int meas_interval = 1;
  std::string outfile;

  if( argc > 1 && argv[1][0] != '-' ){
    config=argv[1];
    interval = stoi(argv[2]);
    nsteps = stoi(argv[3]);
    meas_interval = stoi(argv[4]);
    outfile=argv[5];
    step_size=double(interval)/nsteps;
    std::cout << GridLogMessage << "Loading configuration from " << config << std::endl;
    FieldMetaData header;
    NerscIO::readConfiguration(U, header, argv[1]);
  }
  else {
    SU<Nc>::ColdConfiguration(U);
    config="ColdConfig";
    outfile=config+".h5";
  }

  // RealD plaq = WilsonLoops<Impl>::avgPlaquette(U);


  // wflow<PeriodicGimplD> wflowobs(nsteps, step_size, meas_interval);
  // wflowobs( U );
  LatticeGaugeField Usmear(U);
  
    std::vector<FlowResult> flow_result;
    flow_result.resize(8);

    flow_result[0].description = "Flow time";
    flow_result[1].description = "Plaquette energy density";
    flow_result[2].description = "Clover energy density";
    flow_result[3].description = "Topological charge";
    flow_result[4].description = "Avg Plaquette";
    flow_result[5].description = "Avg Rectangle";
    flow_result[6].description = "Real part of Polyakov loop";
    flow_result[7].description = "Imaginary part of Polyakov loop";

    // put t=0 data (before the flow)
    flow_result[0].data.push_back(0);
    flow_result[1].data.push_back(0);
    flow_result[2].data.push_back(0);
    flow_result[3].data.push_back(WilsonLoops<PeriodicGimplR>::TopologicalCharge(U));
    flow_result[4].data.push_back(WilsonLoops<PeriodicGimplR>::avgPlaquette(U));
    flow_result[5].data.push_back(WilsonLoops<PeriodicGimplR>::avgRectangle(U));
    flow_result[6].data.push_back(WilsonLoops<PeriodicGimplR>::avgPolyakovLoop(U).real());
    flow_result[7].data.push_back(WilsonLoops<PeriodicGimplR>::avgPolyakovLoop(U).imag());

    // WilsonFlow<PeriodicGimplR> wflow(par().epsilon, par().Nstep, 1);
    WilsonFlow<PeriodicGimplR> wflow(step_size, nsteps, meas_interval);
    wflow.resetActions();

    // wflow.addMeasurement(1, [&wflow, &flow_result](int step, RealD t, const LatticeGaugeField &U){
    wflow.addMeasurement(meas_interval, [&wflow, &flow_result](int step, RealD t, const LatticeGaugeField &U){
        flow_result[0].data.push_back(t);
        flow_result[1].data.push_back(wflow.energyDensityPlaquette(t,U));
        flow_result[2].data.push_back(wflow.energyDensityCloverleaf(t,U));
        flow_result[3].data.push_back(WilsonLoops<PeriodicGimplR>::TopologicalCharge(U));
        flow_result[4].data.push_back(WilsonLoops<PeriodicGimplR>::avgPlaquette(U));
        flow_result[5].data.push_back(WilsonLoops<PeriodicGimplR>::avgRectangle(U));
        flow_result[6].data.push_back(WilsonLoops<PeriodicGimplR>::avgPolyakovLoop(U).real());
        flow_result[7].data.push_back(WilsonLoops<PeriodicGimplR>::avgPolyakovLoop(U).imag());

        std::cout << GridLogMessage << "Step " << step << std::endl;
        for (unsigned int i = 0; i < 8; ++i)
        {
            std::cout << GridLogMessage << std::setw(25) << flow_result[i].description << " " << flow_result[i].data.back() << std::endl;
        }
        std::cout << GridLogMessage << std::endl;
    });

    wflow.smear(Usmear, U);

    for (unsigned int i = 1; i < 8; ++i)
    {
      std::cout << GridLogMessage << std::setw(25) << flow_result[i].description << " : " << flow_result[i].data.back() << std::endl;
    }

    Hdf5Writer WR(outfile);
    write(WR,"FlowObservables",flow_result);




  Grid_finalize();
}


