import os
import time
def gen_prm_file(eps , abs1 , rel1 , abs2 , rel2):
    file_prm = open("prm_test_new.txt" , 'w')
    file_prm.write("-0.6\n")
    file_prm.write( str(abs1)+"\n")
    file_prm.write( str(rel1)+"\n")
    file_prm.write( str(abs2)+"\n")
    file_prm.write( str(rel2)+"\n")
    file_prm.write( str(eps) +"\n")
    file_prm.write("0\n")
    file_prm.write("0.1\n")
    file_prm.write("0\n")
    file_prm.close()

if __name__ == '__main__':
    for log10_rel1 in range(1,5):
        for log10_rel2 in range(1,5):
            gen_prm_file(-2.0 , 10**(-log10_rel1) , 0 ,10**(-log10_rel2) , 0  )
	    print( 10**(-log10_rel1) ,10**(-log10_rel2) )
            os.system("./program -c0 -new prm_test_new.txt Data &")
	    time.sleep(0.5)
	    os.system("./program -c0 -old prm_test_new.txt Data &")
	    time.sleep(0.5)