#created by codex

#bash file to run remaining channels automatically, nue_Fe56 not included
#!/usr/bin/env bash

#!/usr/bin/env bash
set -euo pipefail

cd "$HOME/software/geant4/geant4-work/FDSteel/build"

mkdir -p logs

MARLEY="$HOME/software/marley/FDSteel"

channels=(
    nuebar_Fe56
    nue_e
    nuebar_e
    numu_e
    numubar_e
)

layers=(PrimaryIron OuterIron)

# Check all macro files before starting
for layer in "${layers[@]}"; do
    for channel in "${channels[@]}"; do
        for i in {1..8}; do
            macro="macros/marley_scan100k_${layer}_${channel}/marleyScan100k${i}.mac"
            if [[ ! -f "$macro" ]]; then
                echo "Missing macro: $macro" >&2
                exit 1
            fi
        done
    done
done

for layer in "${layers[@]}"; do

    for channel in "${channels[@]}"; do

        case "$channel" in
            nuebar_Fe56)
                config="FDSteel_CC_nuebar_Fe56.js"
                ;;
            *)
                config="FDSteel_eES_Fe56.js"
                ;;
        esac

        export FDSTEEL_CONFIG="$MARLEY/$config"

        if [[ ! -f "$FDSTEEL_CONFIG" ]]; then
            echo "Missing MARLEY config: $FDSTEEL_CONFIG" >&2
            exit 1
        fi

        echo "Starting: $layer / $channel"
        echo "MARLEY config: $FDSTEEL_CONFIG"

        pids=()

        for i in {1..8}; do
            ./exampleB1 \
                "macros/marley_scan100k_${layer}_${channel}/marleyScan100k${i}.mac" \
                > "logs/${layer}_${channel}_${i}.log" 2>&1 &

            pids+=("$!")
        done

        failed=0

        for pid in "${pids[@]}"; do
            if ! wait "$pid"; then
                failed=1
            fi
        done

        if (( failed )); then
            echo "ERROR: At least one job failed for $layer / $channel"
            exit 1
        fi

        echo "Finished: $layer / $channel"

    done

done

echo "All remaining iron simulations completed!"

