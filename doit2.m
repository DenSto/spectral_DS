clear all;
hm_l = 1:0.5:3.5;
dm_l = 1:0.5:5;
for hm = hm_l
	for dm = dm_l
    	spectral_DS2d(hm,dm);
    	clear spectral_DS2d
        clc
	end
end
    
