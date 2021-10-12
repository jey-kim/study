#!/bin/bash

for i in {39..45}; do

if [[ $i -eq 41 ]]
then
  echo "i is 41."
  k=0.18
  echo $k
elif [[ $i -eq 42 ]]
then
  echo "i is 42."
  k=0.19
  echo $k
elif [[ $i -eq 43 ]]
then
  echo "i is 43."
  k=0.22
  echo $k
else
  echo "i is NOT 41,42,nor 43."
  k=0.2
  echo $k
fi

done
