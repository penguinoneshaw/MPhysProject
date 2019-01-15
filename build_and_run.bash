#!/bin/bash
cd build;
cmake3 --build . --config Debug;
cd ..;
nohup build/ProjectModelling data;
