#!/bin/bash
# take_aaron_4q_ops_claude.sh
#
# Fetch Aaron Meyer's 4q nucleon operator tarball (cubic-group irreps, CG "blanks")
# that he sent via LC give/take, pull it here, and extract it locally.
#   1. ssh to LC, `take meyer54` into REMOTE_DIR (on lustre)
#   2. rsync the tarball -> LOCAL_DIR
#   3. list + extract locally (~25MB -> ~700MB)
# All output is tee'd to take_aaron_4q_ops_claude.log.
#
# Usage:
#   ./take_aaron_4q_ops_claude.sh
#   GIVER=otheruser ./take_aaron_4q_ops_claude.sh     # if Aaron's LC username differs
#   SKIP_TAKE=1 ./take_aaron_4q_ops_claude.sh         # file already taken, just pull + extract

set -u

LOG=/mnt/baracuda_14/grid_claude/take_aaron_4q_ops_claude.log
exec > >(tee "${LOG}") 2>&1

REMOTE=matsumoto5@oslic.llnl.gov
GIVER="${GIVER:-meyer54}"
REMOTE_DIR="${REMOTE_DIR:-/p/lustre5/matsumoto5/aaron_4q_ops}"
LOCAL_DIR="${LOCAL_DIR:-/mnt/baracuda_14/grid_claude/baryon_group_theory/aaron_4q_ops}"

# SSH connection multiplexing: authenticate ONCE, reuse the master for every step.
SSH_OPTS="-o ControlMaster=auto -o ControlPath=${HOME}/.ssh/cm-take-%C -o ControlPersist=5m"

echo "Opening master SSH connection to ${REMOTE} (authenticate once) ..."
ssh ${SSH_OPTS} "${REMOTE}" true || { echo "ERROR: ssh failed" ; exit 1 ; }

echo "==== [1] take from ${GIVER} into ${REMOTE_DIR} ===="
if [ "${SKIP_TAKE:-0}" = "1" ]; then
    echo "SKIP_TAKE=1: not running take"
else
    # -t: take may prompt (e.g. overwrite confirmation)
    ssh -t ${SSH_OPTS} "${REMOTE}" "mkdir -p ${REMOTE_DIR} && cd ${REMOTE_DIR} && take ${GIVER}"
    if [ $? -ne 0 ]; then
        echo "ERROR: take failed (wrong giver name? nothing pending?)"
        exit 1
    fi
fi
ssh ${SSH_OPTS} "${REMOTE}" "ls -la ${REMOTE_DIR}"

echo "==== [2] rsync ${REMOTE}:${REMOTE_DIR}/ -> ${LOCAL_DIR}/ ===="
mkdir -p "${LOCAL_DIR}"
rsync -av --progress -e "ssh ${SSH_OPTS}" "${REMOTE}:${REMOTE_DIR}/" "${LOCAL_DIR}/"
if [ $? -ne 0 ]; then
    echo "ERROR: rsync failed"
    exit 1
fi
ssh ${SSH_OPTS} -O exit "${REMOTE}" 2>/dev/null || true

echo "==== [3] extract locally in ${LOCAL_DIR} ===="
cd "${LOCAL_DIR}" || exit 1
ls -la
for tb in *.tar *.tar.gz *.tgz *.tar.bz2 *.tar.xz; do
    [ -e "${tb}" ] || continue
    echo "---- ${tb}: top-level entries ----"
    tar tf "${tb}" | awk -F/ '{print $1"/"$2}' | sort -u | head -50
    echo "---- extracting ${tb} ----"
    tar xf "${tb}"
    if [ $? -ne 0 ]; then
        echo "ERROR: extract of ${tb} failed"
        exit 1
    fi
done

echo "==== summary ===="
du -sh "${LOCAL_DIR}"
find "${LOCAL_DIR}" -maxdepth 4 -type d | sort | head -80
echo "op file count:"
find "${LOCAL_DIR}" -name '*.op' | wc -l
echo "DONE"
