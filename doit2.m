clear all;
hm_l = 1:0.5:3.5;
dm_l = 1:0.5:5;
x=length(hm_l);
y=length(dm_l);
arr=zeros(x*y,3);
arr2=0;
l=1;
for hm = hm_l
	for dm = dm_l
        arr(l,1) = hm;
        arr(l,2) = dm;
        arr(l,3) = dm^2 + hm^2;
    	%spectral_DS2d(hm,dm);
    	%clear spectral_DS2d
        %clc
        l=l+1;
    end
   arr2=sortrows(arr,3)
end
for i = 1:(x*y)
   spectral_DS2d(arr2(i,1),arr2(i,2));
   clear spectral_DS2d
   clc
      
end
    
