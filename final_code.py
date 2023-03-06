import re
import sys
import argparse
import json
import matplotlib.pyplot as plt
import numpy as np
import matplotlib


def parsed_data(data):

    x_data = float(data[2])
    y_data = float(data[3])
    z_data = float(data[4])

    #x_cord
    if x_data > CELL_SIZE/2:
        x_data = x_data - CELL_SIZE

    #y_cord
    if y_data > CELL_SIZE/2:
        y_data = y_data - CELL_SIZE

    #z_cord
    if z_data > CELL_SIZE/2:
        z_data = z_data - CELL_SIZE



    return {
        "atom_id" : data[0],
        "atom_class" : data[1],
        "atom_coordinate" : [x_data,y_data,z_data] 
    }


def get_points(path):
    with open(path, 'r') as myfile:
        final_data = json.load(myfile)

    # pts = [[pt['atom_coordinate'] for pt in pts] for pts in final_data]
    pts = [[{'coordinate':pt['atom_coordinate'],'id' : pt['atom_id'], 'atom_class' : pt['atom_class']} for pt in pts] for pts in final_data]
    return pts


def countX(lst, x):
    count = 0
    for ele in lst:
        if (ele == x):
            count = count + 1
    return count



def get_atom_distribution(data):
	my_distro = {}
	for i in range(MAX_ATOMS_PER_BUCKET+1):
	    my_distro[i] = 0

	# #Only using 1 data
	# all_data = [all_data[0]]
	buckets = {}
	for i in range(NO_OF_BUCKET):
	    buckets[i] = 0
	for each in data:
		datum = each['coordinate']
		x_cord, y_cord, z_cord = datum[0], datum[1], datum[2]
		x_offset = int((x_cord - LEFT_BOUNDARY)/BUCKET_LENGTH)
		y_offset = int((y_cord - LEFT_BOUNDARY)/BUCKET_LENGTH) 
		z_offset = int((z_cord - LEFT_BOUNDARY)/BUCKET_LENGTH)

		z_multiplier = round(NO_OF_BUCKET ** (2/3),2)
		y_multiplier = round(NO_OF_BUCKET ** (1/3),2)

		bucket_no = int((z_multiplier * z_offset + y_multiplier * y_offset + x_offset))
		try:
		    a = buckets[bucket_no]
		    a = a + 1
		    buckets[bucket_no] = a
		except Exception as e:
		    continue

	list_num_per_bucket = []
	for x,y in buckets.items():
	    list_num_per_bucket.append(y)

	unique_elements = list(set(list_num_per_bucket))

	for element in unique_elements:
	    count = countX(list_num_per_bucket,element)
	    my_distro[element] = my_distro[element] + count


	final_frequency = {}
	for key,value in my_distro.items():
	    final_frequency[key] = round(value,2)

	val = []
	for x,y in final_frequency.items():
	    val = val + int(round(y,2)) * [x]


	final_list = []
	for i in val:
	    final_list.append(int(round(i,2)))

	set_list = set(final_list)
	length = len(set_list)


	#Ploting
	from datetime import datetime
	now = datetime.now()
	timestamp = int(now.timestamp())

	try:
	    n, bins, patches = plt.hist(final_list, length, facecolor='blue', alpha=0.5,edgecolor="red",align='mid' )
	    plt.xlabel('Number of Fe atoms in a bin')
	    plt.ylabel('Frequency')
	    plt.title('Distribution of Fe atoms in a cell')
	    plt.savefig('{}_{}'.format(args.output,timestamp))
	except Exception as e:
	    print(e)

#start = mimimum no of atoms per bucket, and end = maximum per bucket
def get_filtered_data(data,start,end):
	start = int(start)
	end = int(end)
	buckets = {}
	for i in range(NO_OF_BUCKET):
	    buckets[i] = 0
	for each in data:
		datum = each['coordinate']
		x_cord, y_cord, z_cord = datum[0], datum[1], datum[2]
		x_offset = int((x_cord - LEFT_BOUNDARY)/BUCKET_LENGTH)
		y_offset = int((y_cord - LEFT_BOUNDARY)/BUCKET_LENGTH) 
		z_offset = int((z_cord - LEFT_BOUNDARY)/BUCKET_LENGTH)

		z_multiplier = round(NO_OF_BUCKET ** (2/3),2)
		y_multiplier = round(NO_OF_BUCKET ** (1/3),2)

		bucket_no = int((z_multiplier * z_offset + y_multiplier * y_offset + x_offset))
		try:
		    a = buckets[bucket_no]
		    a = a + 1
		    buckets[bucket_no] = a
		except Exception as e:
		    continue

	list_num_per_bucket = []
	"""later"""
	return_buckets = []

	for x,y in buckets.items():
	    list_num_per_bucket.append(y)
	    if y >=start and y<=end:
	        return_buckets.append(x)

	out_file = open("filtered_bucket.json", "w")
	  
	json.dump(return_buckets, out_file)
	out_file.close()
	filtered_data = []
	for obj in data:
		datum = obj['coordinate']
		datum_id = obj['id']
		x_cord, y_cord, z_cord = datum[0], datum[1], datum[2]
		x_offset = int((x_cord - LEFT_BOUNDARY)/BUCKET_LENGTH)
		y_offset = int((y_cord - LEFT_BOUNDARY)/BUCKET_LENGTH) 
		z_offset = int((z_cord - LEFT_BOUNDARY)/BUCKET_LENGTH)

		z_multiplier = round(NO_OF_BUCKET ** (2/3),2)
		y_multiplier = round(NO_OF_BUCKET ** (1/3),2)

		bucket_no = int((z_multiplier * z_offset + y_multiplier * y_offset + x_offset))
		if bucket_no in return_buckets:
			# if x_cord > 0:
			# 	x_cord = x_cord - round(CELL_SIZE/2,2)
			# else:
			# 	x_cord = x_cord + round(CELL_SIZE/2,2)
			if y_cord > 0:
				y_cord = y_cord - round(CELL_SIZE/2,2)
			else:
				y_cord = y_cord + round(CELL_SIZE/2,2)
			# if z_cord > 0:
			# 	z_cord = z_cord - round(CELL_SIZE/2,2)
			# else:
			# 	z_cord = z_cord + round(CELL_SIZE/2,2)
			filtered_data.append(
			        {
			            'id' : datum_id,
			            'atom_coordinate' : [round(x_cord,2),round(y_cord,2),round(z_cord,2)] 
			        }
			    )

	if args.input == 'fe.json':
		out_file = open("filtered_fe.json", "w")
	elif args.input == 'mg.json':
		out_file = open("filtered_mg.json", "w")
	elif args.input == 'o.json':
		out_file = open("filtered_o.json", "w")
	else:
		out_file = open("filtered_si.json", "w")

	
	json.dump(filtered_data, out_file)
	out_file.close()


	# f = open("out_filtered.dump", "w")
	# f.write("ITEM: TIMESTEP\n")
	# f.write("0\n")
	# f.write("ITEM: NUMBER OF ATOMS\n")
	# f.write("{}\n".format(len(filtered_data)))
	# f.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
	# f.write("0.0000000000000000e+00 6.8000000000000000e+01 0.0000000000000000e+00\n")
	# f.write("0.0000000000000000e+00 6.8000000000000000e+01 0.0000000000000000e+00\n")
	# f.write("0.0000000000000000e+00 6.8000000000000000e+01 0.0000000000000000e+00\n")
	# f.write("ITEM: ATOMS id type x y z\n")



	# for each in filtered_data:
	#     f.write("{} {} {} {} {}".format(each['id'],"1",str(each['atom_coordinate'][0]),str(each['atom_coordinate'][1]),str(each['atom_coordinate'][2])))
	#     f.write("\n")
	# f.close()

def get_filtered_data_all(data,start,end,atom_class):
	start = int(start)
	end = int(end)
	buckets = {}
	for i in range(NO_OF_BUCKET):
	    buckets[i] = 0
	for each in data:
		print("each=",each)
		if each['atom_class'] != atom_class:
			continue
		datum = each['coordinate']
		x_cord, y_cord, z_cord = datum[0], datum[1], datum[2]
		x_offset = int((x_cord - LEFT_BOUNDARY)/BUCKET_LENGTH)
		y_offset = int((y_cord - LEFT_BOUNDARY)/BUCKET_LENGTH) 
		z_offset = int((z_cord - LEFT_BOUNDARY)/BUCKET_LENGTH)

		z_multiplier = round(NO_OF_BUCKET ** (2/3),2)
		y_multiplier = round(NO_OF_BUCKET ** (1/3),2)

		bucket_no = int((z_multiplier * z_offset + y_multiplier * y_offset + x_offset))
		try:
		    a = buckets[bucket_no]
		    a = a + 1
		    buckets[bucket_no] = a
		except Exception as e:
		    continue

	list_num_per_bucket = []
	"""later"""
	return_buckets = []

	for x,y in buckets.items():
	    list_num_per_bucket.append(y)
	    if y >=start and y<=end:
	        return_buckets.append(x)

	out_file = open("filtered_bucket.json", "w")
	  
	json.dump(return_buckets, out_file)
	out_file.close()
	filtered_data = []
	for obj in data:
		datum = obj['coordinate']
		datum_id = obj['id']
		x_cord, y_cord, z_cord = datum[0], datum[1], datum[2]
		x_offset = int((x_cord - LEFT_BOUNDARY)/BUCKET_LENGTH)
		y_offset = int((y_cord - LEFT_BOUNDARY)/BUCKET_LENGTH) 
		z_offset = int((z_cord - LEFT_BOUNDARY)/BUCKET_LENGTH)

		z_multiplier = round(NO_OF_BUCKET ** (2/3),2)
		y_multiplier = round(NO_OF_BUCKET ** (1/3),2)

		bucket_no = int((z_multiplier * z_offset + y_multiplier * y_offset + x_offset))
		if bucket_no in return_buckets:
			# if x_cord > 0:
			# 	x_cord = x_cord - round(CELL_SIZE/2,2)
			# else:
			# 	x_cord = x_cord + round(CELL_SIZE/2,2)
			if y_cord > 0:
				y_cord = y_cord - round(CELL_SIZE/2,2)
			else:
				y_cord = y_cord + round(CELL_SIZE/2,2)
			# if z_cord > 0:
			# 	z_cord = z_cord - round(CELL_SIZE/2,2)
			# else:
			# 	z_cord = z_cord + round(CELL_SIZE/2,2)
			filtered_data.append(
			        {
			            'id' : datum_id,
			            'atom_coordinate' : [round(x_cord,2),round(y_cord,2),round(z_cord,2)],
			            'atom_class' : obj['atom_class'] 
			        }
			    )

	if atom_class == "1":
		out_file = open("cluster_within_fe.json", "w")
	elif atom_class == "2":
		out_file = open("cluster_within_mg.json", "w")
	elif atom_class == "3":
		out_file = open("cluster_within_o.json", "w")
	else:
		out_file = open("cluster_within_si.json", "w")

	
	json.dump(filtered_data, out_file)
	out_file.close()

	f_out = open("merged_cluster.dump", "w")
	f_out.write("ITEM: TIMESTEP\n")
	f_out.write("0\n")
	f_out.write("ITEM: NUMBER OF ATOMS\n")
	f_out.write("{}\n".format(str(len(filtered_data))))
	f_out.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
	f_out.write("0.0000000000000000e+00 6.8000000000000000e+01 0.0000000000000000e+00\n")
	f_out.write("0.0000000000000000e+00 6.8000000000000000e+01 0.0000000000000000e+00\n")
	f_out.write("0.0000000000000000e+00 6.8000000000000000e+01 0.0000000000000000e+00\n")
	f_out.write("ITEM: ATOMS id type x y z\n")


	for each in filtered_data:
		f_out.write("{} {} {} {} {}".format(each['id'],each['atom_class'],str(each['atom_coordinate'][0]),str(each['atom_coordinate'][1]),str(each['atom_coordinate'][2])))
		f_out.write("\n")

	f_out.close()




def merge_data(merge_val):
	total = 0 
	for i in merge_val:
		if i == "1":
			f = open("filtered_fe.json","r")
		elif i == "2":
			f = open("filtered_mg.json","r")
		elif i == "3":
			f = open("filtered_o.json","r")
		else:
			f = open("filtered_si.json","r")
		filtered_data = json.load(f)
		total = total + len(filtered_data)

	f_out = open("out_filtered.dump", "w")
	f_out.write("ITEM: TIMESTEP\n")
	f_out.write("0\n")
	f_out.write("ITEM: NUMBER OF ATOMS\n")
	f_out.write("{}\n".format(str(total)))
	f_out.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
	f_out.write("0.0000000000000000e+00 6.8000000000000000e+01 0.0000000000000000e+00\n")
	f_out.write("0.0000000000000000e+00 6.8000000000000000e+01 0.0000000000000000e+00\n")
	f_out.write("0.0000000000000000e+00 6.8000000000000000e+01 0.0000000000000000e+00\n")
	f_out.write("ITEM: ATOMS id type x y z\n")
	for i in merge_val:
		if i == "1":
			f = open("filtered_fe.json","r")
		elif i == "2":
			f = open("filtered_mg.json","r")
		elif i == "3":
			f = open("filtered_o.json","r")
		else:
			f = open("filtered_si.json","r")
		filtered_data = json.load(f)

		for each in filtered_data:
			f_out.write("{} {} {} {} {}".format(each['id'],str(i),str(each['atom_coordinate'][0]),str(each['atom_coordinate'][1]),str(each['atom_coordinate'][2])))
			f_out.write("\n")
	f_out.close()

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Parser')
	parser.add_argument('-n','--no_of_atoms', help='No of atoms',required=True)
	parser.add_argument('-s','--cell_size', help='Maximimum cell size',required=True)
	parser.add_argument('-b','--no_of_buckets', help='No of buckets',required=True)
	parser.add_argument('-p','--percent', help='nth simulation',required=True)
	parser.add_argument('-i','--input', help='Input Json File',required=False)
	parser.add_argument('-o','--output', help='Output File Name',required=True)
	parser.add_argument('-k','--to_do', help='What to perform',required=True)
	parser.add_argument('-f','--filter', help='Filter Values',required=False)
	parser.add_argument('-m','--merge', help='What data to merge',required=False)
	args = parser.parse_args()


	NO_OF_ATOMS = int(args.no_of_atoms)
	CELL_SIZE = int(args.cell_size)
	NO_OF_BUCKET = int(args.no_of_buckets)

	BOUNDARY = (round(-1 * CELL_SIZE/2,2),round(CELL_SIZE/2,2))
	LEFT_BOUNDARY = BOUNDARY[0]
	RIGHT_BOUNDARY = BOUNDARY[1]

	BUCKET_LENGTH = CELL_SIZE/ (NO_OF_BUCKET**(1/3)) + 0.1 #adding val for boundary conditions

	AVG_ATOMS_PER_BUCKET = int(NO_OF_ATOMS/NO_OF_BUCKET)
	#setting max atom a bucket can have for x-axis
	MAX_ATOMS_PER_BUCKET = 2 * AVG_ATOMS_PER_BUCKET

	if args.input:
		all_data = get_points(args.input)
	else:
		all_data = []
	#%percentile of data
	n_data =int(args.percent)

	if n_data > 1 or n_data <0:
		val = len(all_data) - 1
	else:
		#Percentage of nth data
		val = int(n_data * len(all_data)) - 1

	if args.to_do == "distribution":
		get_atom_distribution(all_data[val])
	elif args.to_do == "filter":
		f = args.filter
		if not f:
			print("Filter Parameter not passed")
			exit(1)
		start = f.split('-')[0]
		end = f.split('-')[1]
		get_filtered_data(all_data[val],start,end)
	elif args.to_do == "merge":
		merge_val = args.merge.split("-")
		merge_data(merge_val)
	elif args.to_do == "nested_cluster":
		f = args.filter
		if not f:
			print("Filter Parameter not passed")
			exit(1)
		start = f.split('-')[0]
		end = f.split('-')[1]
		atom_class = f.split('-')[2]
		get_filtered_data_all(all_data[val],start,end,atom_class)
	else:
		exit(1)





























