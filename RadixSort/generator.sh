#!/bin/bash
# uniform_max=(10 100 1000 5000 7000 8000 10000 15000 20000 50000 100000 1000000 10000000 100000000 1000000000)
# testcase=(4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 37 38 39 40 41 42 43)
# testcase=(30 31 32 33 34 35)
testcase=(21 22 23 24 25 26 27)
for ((i=0;i<${#testcase[@]};i++))
do
    nohup numactl -i all ./radixSort ${testcase[i]} >> radixSort_8_8pair_192.txt
    # echo "Uniform Distribution ${i}" >> test_parlay_sample_0915.txt
    # echo "Key range  = (0, ${uniform_max[i]})" >> test_parlay_sample_0915.txt
    # nohup numactl -i all ./test_parlay ${NUM_ELEMENTS} uniform ${uniform_max[i]} >> test_parlay_sample_0915.txt
done

