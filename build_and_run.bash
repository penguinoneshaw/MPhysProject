#!/bin/bash
cd build;
cmake3 --build . --config Release || exit 1;
cd ..;
stat output && mv output output_prev;
OMP_NUM_THREADS=4 build/ProjectModelling data || exit 1;
tar czf output.tar.gz output
echo "Results from data run done `date`" | mail -a output.tar.gz -s "project run complete" s1410767
exit 0;
