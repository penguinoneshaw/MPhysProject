#!/bin/bash
cd build;
cmake3 --build . --config Debug;
cd ..;
mv output output_prev
nohup build/ProjectModelling data;
