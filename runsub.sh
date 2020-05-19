t=gmask3
d0=.5
d1=5
r=False
#{7..14}
for i in 11
do
	python subcubes.py -fitcat EAGLE -dmin $d0 -dmax $d1 -snap $i -type $t -remove $r& 
done

