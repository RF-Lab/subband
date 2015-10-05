clc, clear all;

curPath = pwd() ;
cd('.\tsim') ;
modelPath = pwd() ;
cd( curPath ) ;
addpath(modelPath) ;

fs = 5;
fd = 20 ;
N = 10;

freq_range = 1; w1 = (fs + freq_range)/fd * 2 *pi; 

F0_order2 = get_sb_matrix_2(N, w1, 0); 
F1_order2 = get_sb_matrix_2(N, w1, 1); 
F2_order2 = get_sb_matrix_2(N, w1, 2); 

s = cos(2*pi*fs/fd * (0:N-1));

rxx2_sb2 = get_acf_from_sb_matrix_2(N, F0_order2, F1_order2, F2_order2, s);                
sub_pi = ar_model(rxx2_sb2) ;
[poles1, omega0_sub_2nd, Hjw0_1_pi] = get_ar_pole(sub_pi); 
fs_sub_2nd = omega0_sub_2nd*fd/2/pi;

assert((fs - fs_sub_2nd) < 0.1, 'Assert. Subband estimation doesnt work');

% remove model path
rmpath(modelPath) ;