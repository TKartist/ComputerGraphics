#!/bin/bash
mkdir results &> /dev/null

rm -r ./results/* &> /dev/null

g++ main.cpp -o main.out

for i in {00..120}
do
./main.out "./results/result_${i}.ppm" ${i} && ffmpeg -y -i "./results/result_${i}.ppm" "./results/result_${i}.jpg"
done

ffmpeg -y -framerate 60 -pattern_type glob -i './results/*.jpg' -c:v libx264 -pix_fmt yuv420p ./result.mp4

ffmpeg -y -framerate 60 -pattern_type glob -i './results/*.jpg' -loop 0 ./result.gif
