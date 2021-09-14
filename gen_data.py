def proc_data_error(name_file_error):
    file = open(r'OUT_Error')
    data = [  ]
    set_s = set()
    for line in file:
        substrs = line.split()
        x = float(substrs[0])
        s = float(substrs[1])
        t = float(substrs[2])
        e = float(substrs[3])
        arg = float(substrs[4])
        M = float(substrs[5])
        R = float(substrs[6])
        tmp =[]
        tmp.append(s)
        tmp.append(x)
        tmp.append(arg)
        tmp.append(R)
        data.append(tmp)
        set_s.add(s)
    for i in set_s:
	if abs( i + 0.6) < 0.01:
	        files = open("data_s_" + str(i) ,'w')
	        for arr in data:
	        	if arr[0] == i :
				for j in range(3):
					files.write(str(arr[1+j]))
					files.write('\t')
				files.write('\n')
    data.sort()
    load_file = open("new_data" , 'w')
    for line in data:
        for i in line:
            load_file.write(str(i))
            load_file.write('\t')
        load_file.write('\n')

def proc_data_OUT(name_file_OUT):
	file_OUT = open(name_file_OUT)
	data_rho = []
	for line in file_OUT:
		arr_data = line.split()
		tmp_data_rho = []
		for count in range(4):
			tmp_data_rho.append(float(arr_data.pop(0))) # \epsilon \rho lower_bounds upp_bounds
		data_rho.append(tmp_data_rho)
		file_data_for_first_int = open("data_for_eps_" + str(tmp_data_rho[0]) , 'w' )
		tmp_arr = []
		for i in arr_data:
			tmp_arr.append(i)
			if len(tmp_arr) == 4:
				for word in tmp_arr:
					file_data_for_first_int.write(word + '\t')
				file_data_for_first_int.write('\n')
				tmp_arr = []		
	file_data_rho = open("data_rho" , 'w')
	for line in data_rho:
		for word in line:
			file_data_rho.write(str(word) + '\t')
		file_data_rho.write('\n')			
				

if __name__ == '__main__':
	proc_data_OUT("OUT_new4")
    #fig , ax = plt.subplots()
    #plt.plot(new_data[0] , new_data[1] , label = "arg" )
    #plt.plot(new_data[0] , new_data[3] , label = "R" )
    #plt.legend()
    #plt.grid()
    #plt.show()
    #fig.save("graph.png")
