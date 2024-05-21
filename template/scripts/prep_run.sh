#!/bin/bash


# Set the source code directory
export JULIADEPOT="/work/n02/n02/aleey/wavi/WAVIhpc/.julia_AlexEKI"


# Empty the run directory - but first make sure it exists!
if [ -d "../run" ]; then
  cd ../run
  rm -rf *
else
  echo 'There is no run directory'
  exit 1
fi

# Link everything from the input directory
#ln -s ../input/* .

# Copy everything from the input directory
cp -f ../input/* .

# Deep copy of the driver file
#rm -f driver.jl
#§cp -f ../input/driver.jl .

# Deep copy of the current version of wavi source directory
mkdir src
cp -rf $JDEPOT/packages/WAVI/. src/.

# Deep copy of any pickups (so they don't get overwritten in input/)
rm -f pickup*
cp -f ../input/pickup* . 2>/dev/null

# Deep copy of any functions in the code directory
mkdir code
cp -f ../code/* code/.

# Soft link the driver code (only necessary if you're pulling in from multiple scripts!)
#ln -s ../../../utilities/sub_driver.jl .
#cp -s ../../../utilities/sub_driver.jl .

#Link the image container
ln -s $IMGPATH .
#cp -rf $IMGPATH .

#change to run directory
#cd ../run
