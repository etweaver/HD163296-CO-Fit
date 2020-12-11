#actually run radmc in its folder to generate the temperature
#this is called when an ensemble generation is completed
#it should only ever be called internally
cd radmc
radmc3d mctherm > radmc3dOutput.log
./extract_temp_3d dust_temperature.dat junk.txt
cd ..
#obviously using andrea's program to rearrange the radmc output is wasteful
#and unneccesary. I'll remove it eventually
