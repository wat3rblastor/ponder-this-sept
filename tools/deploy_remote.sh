#!/bin/bash
# Drive a rented multi-GPU machine from this box.
#
#   tools/deploy_remote.sh setup HOST PORT          copy the repo, build, benchmark
#   tools/deploy_remote.sh start HOST PORT [RATIO]  replan, split RATIO:1 remote:local, launch both
#   tools/deploy_remote.sh pull  HOST PORT          loop: fetch remote logs every 30 s
#   tools/deploy_remote.sh stop  HOST PORT          stop the remote engines (checkpointed)
#
# RATIO = remote speed / this box's speed (default 40). Remote results land in
# experiments/remote-HOST/, where tools/autopromote.sh picks up any new record
# and tools/plan_units.py counts the units as done.
cd "$(dirname "$0")/.." || exit 1
MODE=$1; HOST=$2; PORT=$3; RATIO=${4:-40}
SSHO="-p $PORT -o StrictHostKeyChecking=accept-new -o BatchMode=yes -o ConnectTimeout=20"
R="ssh $SSHO root@$HOST"
RD=/root/ponder-this-sept
E=experiments/2026-10-03-cuda
case $MODE in
setup)
  tar czf - --exclude=build --exclude=.git --exclude='experiments/*/*.log' . | $R "mkdir -p $RD && cd $RD && tar xzf -" || exit 1
  $R "cd $RD && export PATH=\$PATH:/usr/local/cuda/bin && nproc && tools/build_here.sh"
  ;;
start)
  # stop the local engine at a unit boundary so its finished units are on disk
  pkill -x run_plan.sh; pkill -INT -x apsearch_cuda; while pgrep -x apsearch_cuda >/dev/null; do sleep 1; done
  STAMP=$(date +%H%M%S)
  MODCAP=2e16 python3 tools/plan_units.py --budget-res ${BUDGET:-2e16} --kmax 600000 --smax 3000 \
      --done-glob "experiments/*/*.jsonl" --out $E/plan_$STAMP.txt || exit 1
  awk -v r=$((RATIO + 1)) 'NR % r == 0' $E/plan_$STAMP.txt > $E/plan_${STAMP}_local.txt
  awk -v r=$((RATIO + 1)) 'NR % r != 0' $E/plan_$STAMP.txt > $E/plan_${STAMP}_remote.txt
  scp -q -P $PORT -o BatchMode=yes $E/plan_${STAMP}_remote.txt root@$HOST:$RD/plan_remote.txt || exit 1
  $R "cd $RD && tools/multi_gpu.sh plan_remote.txt r$STAMP"
  setsid nohup tools/run_plan.sh $E/plan_${STAMP}_local.txt $E/l$STAMP.jsonl 58 > /dev/null 2>&1 < /dev/null &
  echo "local engine restarted on its share; now run: tools/deploy_remote.sh pull $HOST $PORT"
  ;;
pull)
  L=experiments/remote-$HOST; mkdir -p $L
  while true; do
    $R "cd $RD/experiments/remote && tar czf - . 2>/dev/null" | tar xzf - -C $L 2>/dev/null
    sleep 30
  done
  ;;
stop)
  $R "pkill -INT -x apsearch_cuda; sleep 5; pgrep -x apsearch_cuda | wc -l"
  ;;
*) echo "usage: $0 setup|start|pull|stop HOST PORT [RATIO]"; exit 2;;
esac
