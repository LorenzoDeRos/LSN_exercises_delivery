CC = g++

MolDyn_NVE.exe : MolDyn_NVE.cpp
	$(CC) -o MolDyn_NVE.exe MolDyn_NVE.cpp
clean :
	rm output_epot.out
	rm output_ekin.out
	rm output_etot.out
	rm output_temp.out
	rm ave_epot.out
	rm ave_ekin.out
	rm ave_etot.out
	rm ave_temp.out
	rm gdir_ave.out
	rm frames/*.xyz
import_liquid :
	cp input.liquid input.dat
	cp config.liquid config.final
	cp config_old.liquid config_old.final
import_solid :
	cp input.solid input.dat
	cp config.solid config.final
	cp config_old.solid config_old.final
import_gas_argon :
	cp input_argon.gas input.dat
	cp config_argon.gas config.final
	cp config_argon_old.gas config_old.final
import_liquid_argon :
	cp input_argon.liquid input.dat
	cp config_argon.liquid config.final
	cp config_argon_old.liquid config_old.final
import_solid_argon :
	cp input_argon.solid input.dat
	cp config_argon.solid config.final
	cp config_argon_old.solid config_old.final
save_liquid :
	cp output_epot.out liquid_measures_epot.dat
	cp output_ekin.out liquid_measures_ekin.dat
	cp output_etot.out liquid_measures_etot.dat
	cp output_temp.out liquid_measures_temp.dat
	cp ave_epot.out liquid_ave_epot.dat
	cp ave_ekin.out liquid_ave_ekin.dat
	cp ave_etot.out liquid_ave_etot.dat
	cp ave_temp.out liquid_ave_temp.dat
save_solid :
	cp output_epot.out solid_measures_epot.dat
	cp output_ekin.out solid_measures_ekin.dat
	cp output_etot.out solid_measures_etot.dat
	cp output_temp.out solid_measures_temp.dat
	cp ave_epot.out solid_ave_epot.dat
	cp ave_ekin.out solid_ave_ekin.dat
	cp ave_etot.out solid_ave_etot.dat
	cp ave_temp.out solid_ave_temp.dat
save_gas_argon :
	cp output_epot.out gas_argon_measures_epot.dat
	cp output_ekin.out gas_argon_measures_ekin.dat
	cp output_etot.out gas_argon_measures_etot.dat
	cp output_temp.out gas_argon_measures_temp.dat
	cp ave_epot.out gas_argon_ave_epot.dat
	cp ave_ekin.out gas_argon_ave_ekin.dat
	cp ave_etot.out gas_argon_ave_etot.dat
	cp ave_temp.out gas_argon_ave_temp.dat
save_liquid_argon :
	cp output_epot.out liquid_argon_measures_epot.dat
	cp output_ekin.out liquid_argon_measures_ekin.dat
	cp output_etot.out liquid_argon_measures_etot.dat
	cp output_temp.out liquid_argon_measures_temp.dat
	cp ave_epot.out liquid_argon_ave_epot.dat
	cp ave_ekin.out liquid_argon_ave_ekin.dat
	cp ave_etot.out liquid_argon_ave_etot.dat
	cp ave_temp.out liquid_argon_ave_temp.dat
save_solid_argon :
	cp output_epot.out solid_argon_measures_epot.dat
	cp output_ekin.out solid_argon_measures_ekin.dat
	cp output_etot.out solid_argon_measures_etot.dat
	cp output_temp.out solid_argon_measures_temp.dat
	cp ave_epot.out solid_argon_ave_epot.dat
	cp ave_ekin.out solid_argon_ave_ekin.dat
	cp ave_etot.out solid_argon_ave_etot.dat
	cp ave_temp.out solid_argon_ave_temp.dat