#!/bin/bash

echo "Git start..."
git add .
str=$1
if ! [ -n "$1" ]; then
    str="Update"
fi
git commit -m $str
git push origin master
echo "Git end."