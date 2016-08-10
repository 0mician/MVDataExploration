#!/bin/bash

DIR=$(pwd)

while inotifywait -r -e modify $DIR; do
    makeglossaries report
    bibtex report
    make 
done
