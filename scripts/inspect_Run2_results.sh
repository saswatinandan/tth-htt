#!/bin/bash

declare -A RESULTS
declare -A RESULT

RESULT["0l_2tau"]='/local/karl/ttHAnalysis/2016/2020Apr07'
RESULT["1l_1tau"]='/local/karl/ttHAnalysis/2016/2020Apr07'
RESULT["1l_2tau"]='/local/karl/ttHAnalysis/2016/2020Apr07'
RESULT["2l_2tau"]='/local/karl/ttHAnalysis/2016/2020Mar18'
RESULT["2lss_1tau"]='/local/karl/ttHAnalysis/2016/2020Apr07'
RESULT["2los_1tau"]='/local/karl/ttHAnalysis/2016/2020Mar18'
RESULT["2lss"]='/local/karl/ttHAnalysis/2016/2020Apr09'
RESULT["3l"]='/local/karl/ttHAnalysis/2016/2020Mar19'
RESULT["3lctrl"]='/local/karl/ttHAnalysis/2016/2020Mar19'
RESULT["3l_1tau"]='/local/karl/ttHAnalysis/2016/2020Mar18'
RESULT["4l"]='/local/karl/ttHAnalysis/2016/2020Mar19'
RESULT["4lctrl"]='/local/karl/ttHAnalysis/2016/2020Mar19'
ERA=2016
for CHANNEL in "${!RESULT[@]}"; do
  RESULTS[$ERA,$CHANNEL]=${RESULT[$CHANNEL]}
done

RESULT["0l_2tau"]='/local/karl/ttHAnalysis/2017/2020Apr07'
RESULT["1l_1tau"]='/local/karl/ttHAnalysis/2017/2020Apr07'
RESULT["1l_2tau"]='/local/karl/ttHAnalysis/2017/2020Apr07'
RESULT["2l_2tau"]='/local/karl/ttHAnalysis/2017/2020Mar18'
RESULT["2lss_1tau"]='/local/karl/ttHAnalysis/2017/2020Apr07'
RESULT["2los_1tau"]='/local/karl/ttHAnalysis/2017/2020Mar18'
RESULT["2lss"]='/local/karl/ttHAnalysis/2017/2020Apr09'
RESULT["3l"]='/local/karl/ttHAnalysis/2017/2020Mar19'
RESULT["3lctrl"]='/local/karl/ttHAnalysis/2017/2020Mar19'
RESULT["3l_1tau"]='/local/karl/ttHAnalysis/2017/2020Mar18'
RESULT["4l"]='/local/karl/ttHAnalysis/2017/2020Mar19'
RESULT["4lctrl"]='/local/karl/ttHAnalysis/2017/2020Mar19'
ERA=2017
for CHANNEL in "${!RESULT[@]}"; do
  RESULTS[$ERA,$CHANNEL]=${RESULT[$CHANNEL]}
done

RESULT["0l_2tau"]='/local/karl/ttHAnalysis/2018/2020Apr07'
RESULT["1l_1tau"]='/local/karl/ttHAnalysis/2018/2020Apr07'
RESULT["1l_2tau"]='/local/karl/ttHAnalysis/2018/2020Apr07'
RESULT["2l_2tau"]='/local/karl/ttHAnalysis/2018/2020Mar18'
RESULT["2lss_1tau"]='/local/karl/ttHAnalysis/2018/2020Apr07'
RESULT["2los_1tau"]='/local/karl/ttHAnalysis/2018/2020Mar18'
RESULT["2lss"]='/local/karl/ttHAnalysis/2018/2020Apr09'
RESULT["3l"]='/local/karl/ttHAnalysis/2018/2020Mar19'
RESULT["3lctrl"]='/local/karl/ttHAnalysis/2018/2020Mar19'
RESULT["3l_1tau"]='/local/karl/ttHAnalysis/2018/2020Mar18'
RESULT["4l"]='/local/karl/ttHAnalysis/2018/2020Mar19'
RESULT["4lctrl"]='/local/karl/ttHAnalysis/2018/2020Mar19'
ERA=2018
for CHANNEL in "${!RESULT[@]}"; do
  RESULTS[$ERA,$CHANNEL]=${RESULT[$CHANNEL]}
done

ERAS="2016 2017 2018"
OUTPUT_DIR="."

show_help() {
  THIS_SCRIPT=$0;
  echo "Usage: $(basename $THIS_SCRIPT) [-e <era 1> [<era 2>[ <era 3> ]]] [-o <output directory>]";
  exit 0;
}

while getopts "h?e:o:" opt; do
  case "${opt}" in
  h|\?) show_help
        ;;
  e) ERAS=${OPTARG}
        ;;
  o) OUTPUT_DIR=${OPTARG}
        ;;
  esac
done

if [ ! -d $OUTPUT_DIR ]; then
  echo "No such directory: $OUTPUT_DIR";
  exit 1;
fi

declare -A INSPECT_RLE_ARGS
EVENT_YIELDS_ARG=""
for ERA in $ERAS; do
  INSPECT_RLE_ARGS_ERA="";
  for CHANNEL in 0l_2tau 1l_1tau 1l_2tau 2l_2tau 2lss 2lss_1tau 2los_1tau 3l 3lctrl 3l_1tau 4l 4lctrl; do
    if [ ! -z "${RESULTS[$ERA,$CHANNEL]}" ]; then
      EVENT_YIELDS_ARG="${EVENT_YIELDS_ARG} ${RESULTS[$ERA,$CHANNEL]}/histograms/$CHANNEL";
      INSPECT_RLE_ARGS_ERA="$INSPECT_RLE_ARGS_ERA ${RESULTS[$ERA,$CHANNEL]}/output_rle/$CHANNEL";
    fi
  done
  INSPECT_RLE_ARGS[$ERA]=$INSPECT_RLE_ARGS_ERA
done

YIELDS_ERA=$(echo $ERAS | tr ' ' '_');
YIELDS_OUTPUT="$OUTPUT_DIR/event_yields_ERA.log";
INSPECTION_OUTPUT="$OUTPUT_DIR/rle_inspection_ERA.log";
echo "Finding event yields for era(s): $ERAS"
print_event_yields.py -i $EVENT_YIELDS_ARG &> "${YIELDS_OUTPUT//ERA/$YIELDS_ERA}";

for ERA in $ERAS; do
  echo "Inspecting RLE numbers for era $ERA";
  inspect_rle_numbers.py -i ${INSPECT_RLE_ARGS[$ERA]} -v &> "${INSPECTION_OUTPUT//ERA/$ERA}";
done
