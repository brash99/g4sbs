#!/bin/sh
echo 'Running jobs:'
squeue -u $USER | grep R | grep g4sbs | wc -l
echo 'Pending jobs:'
squeue -u $USER | grep PD | grep g4sbs | wc -l
