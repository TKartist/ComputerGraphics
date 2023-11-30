
g++ part2Large.cpp --std=c++11 -O3

declare -a execution_times

measure_execution_time() {
    local start_time
    local end_time
    local duration

    for i in {00..09}
    do
    start_time=$(date +%s.%N)

    ./a.out

    end_time=$(date +%s.%N)
    execution_times+=($start_time $end_time)
    done

    g++ main.cpp --std=c++11 -O3

    for i in {00..09}
    do
    start_time=$(date +%s.%N)
    ./a.out

    end_time=$(date +%s.%N)
    execution_times+=($start_time $end_time)
    done

    g++ part2.cpp --std=c++11 -O3

    for i in {00..09}
    do
    start_time=$(date +%s.%N)
    ./a.out

    end_time=$(date +%s.%N)
    execution_times+=($start_time $end_time)
    done
}

# Example usage
measure_execution_time
python analysis.py "${execution_times[@]}"

# # Run another command or code block
# measure_execution_time
# echo "Execution Time: ${execution_times[1]} seconds"