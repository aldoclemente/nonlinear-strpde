#!/bin/bash

mkdir -p qsuboutput

qsub -o qsuboutput/output-app-amy-negative.txt -e qsuboutput/error--app-amy-negative.txt -N app-amy-negative run_app-amy-negative.qsub.sh
