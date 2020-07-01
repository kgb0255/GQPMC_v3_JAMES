#!bin/bash
#conda activate gqp

num_walkers=40
sim='lgal'
spec_or_photo=specphoto
noise=bgs0_legacy
model=emulator

python get_data.py /global/cscratch1/sd/kgb0255/gqp_mc/mini_mocha/ispeculator/james/emulator_30x/ \
	$num_walkers $sim $spec_or_photo $noise $model
