#!/bin/bash
#cd /
#cd opt/alps/

aValues=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0)


aStr="aValue"

for ii in "${aValues[@]}" 
	do
		cat tempMf | sed -e "s,$aStr,$ii,g" > Makefile
		cat head/tempPara | sed -e "s,$aStr,$ii,g" > head/Parameters.h
#		rm tmpXYZD tmpXYZDH tmpXYZDHJX tmpXYZDHJXJY xyz
#		cd bin
		make
#		./arata$ii #serial programming
		gnome-terminal -e "bash -c \"./arata$ii; exec bash\"" #open new terminal and run ./arata$ii
#		mpirun -np 22 loop --mpi --write-xml ../tez/xyz.in.xml
#		#sudo ./loop --write-xml ../tez/xyz.in.xml
#		cd /opt/alps/tez 
#		tar -cf tez.d=${d}.h=${h}.Jx=${Jx[$j]}.Jy=${Jy[$j]}.Jz=${Jz[$j]}.tar *.* 
#	        cp -r tez.d=${d}.h=${h}.Jx=${Jx[$j]}.Jy=${Jy[$j]}.Jz=${Jz[$j]}.tar ../results
#		rm -rf *.*						
#		cd ..		
	done

echo "Simulasyonlar tamamlandi, yeni hesaplama koyabilirsiniz." | mail -s "Simulasyonlar tamamlandi." saba.karakas@gmail.com
