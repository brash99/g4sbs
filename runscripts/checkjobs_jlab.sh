#!/bin/sh
echo 'Running jobs:'
squeue -u brash | grep R | grep g4sbs | wc -l
echo 'Pending jobs:'
squeue -u brash | grep PD | grep g4sbs | wc -l
