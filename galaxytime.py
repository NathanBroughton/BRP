import h5py
import numpy as np
import matplotlib.pyplot as plt
import eagleSqlTools as sql
from read_eagle import EagleSnapshot
G = 4.302e-3 #pc*Ms^-1*(km/s)^2

def findr(x1,y1,z1,x2,y2,z2):
	#print x1,x2,y1,y2,z1,z2
	distance = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
	#print distance
	#if(distance > Rmax):
	#	Rmax = distance
	return distance

def findRmax(x1,y1,z1,x2,y2,z2,Rmax):
        #print x1,x2,y1,y2,z1,z2
        distance = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
        #print distance
        if(distance > Rmax):
               Rmax = distance
        return distance


def Vectorlength(x,y,z):
	V = np.sqrt(x**2 + y**2 + z**2)
	return V

def Tdym(R,M, V):
	#V = np.sqrt((G*M)/R)
	Rkm = R * 10**6 * 3.086e13
	t = (2*np.pi*Rkm)/V
	return t

def Tfric(R,r,M,m, V):
	T = Tdym(R, M, V) #sec
	print R,r,M,m,T
	t = (1/(8*np.log(R*10**6/r*10**6))) * (M/m) * T
	t = t / (3600*24*365*10**9) # sec -> Gyr
	#print t
	return t 

def read_header(index):
	if indices[index+1] == "19":
		a = 0.5
        if indices[index+1] == "20":
                a = 0.54
        if indices[index+1] == "21":
                a = 0.58
        if indices[index+1] == "22":
                a = 0.62
        if indices[index+1] == "23":
                a = 0.67
        if indices[index+1] == "24":
                a = 0.73
        if indices[index+1] == "25":
                a = 0.79
	if indices[index+1] == "26":
		a = 0.85
	if indices[index+1] == "27":
       		a = 0.91         
	f = h5py.File('/disks/eagle/L0100N1504/REFERENCE/data/snapshot_027_z000p101/snap_027_z000p101.0.hdf5', 'r')
	h = f['Header'].attrs.get('HubbleParam')
        boxsize = f['Header'].attrs.get('Boxsize')
        f.close()
        return a,h, boxsize

class Findmass:
	def __init__(self, gn, sgn, centre,index,load_region_length=0.5):
		self.a, self.h, self.boxsize = read_header(index)
		self.centre = centre
		self.read_galaxy(gn,sgn,centre,index,load_region_length)
	def read_galaxy(self,gn, sgn, centre,index,load_region_length):
		print indices[index]
		indices2 = "%i %s %s " %(indices[index] + 6,indices[index+1],indices[index+2])
		fil.write(indices2)

		data = {}
		eagle_data = EagleSnapshot('/disks/eagle/L0100N1504/REFERENCE/data/snapshot_027_z000p101/snap_027_z000p101.0.hdf5')
	       	centre *= self.h
		
                region = np.array([(centre[0]-0.5*load_region_length), (centre[0]+0.5*load_region_length), (centre[1]-0.5*load_region_length), (centre[1]+0.5*load_region_length), (centre[2]-0.5*load_region_length),(centre[2]+0.5*load_region_length)])

                eagle_data.select_region(*region)
		             
		f = h5py.File('/disks/eagle/L0100N1504/REFERENCE/data/snapshot_027_z000p101/snap_027_z000p101.0.hdf5', 'r')
                constants = f['Constants']
		Mpc = constants.attrs['CM_PER_MPC']
		totalmass = 0
		Rmax = 0
		r = 0	
		Velocity = 0
                


		tmp = eagle_data.read_dataset(0, 'Coordinates')
		cgs = f['PartType%i/%s'%(0,'Coordinates')].attrs.get('CGSConversionFactor')
        	aexp = f['PartType%i/%s'%(0,'Coordinates')].attrs.get('aexp-scale-exponent')
                hexp = f['PartType%i/%s'%(0,'Coordinates')].attrs.get('h-scale-exponent')

		#print cgs, aexp, hexp
		centre2 = np.multiply(centre, cgs / Mpc * self.a**aexp *self.h**hexp)
		print self.a, self.h, cgs, aexp, hexp
		for j in range(len(myData)):
			Velocity += Vectorlength(myData[j][10], myData[j][11], myData[j][12])
			myData[j][2] = np.multiply(myData[j][2], cgs / Mpc * self.a**aexp *self.h**hexp *self.h)
                 	myData[j][3] = np.multiply(myData[j][3], cgs / Mpc * self.a**aexp *self.h**hexp *self.h)
                       	myData[j][4] = np.multiply(myData[j][4], cgs / Mpc * self.a**aexp *self.h**hexp *self.h)


		for j in range(len(myData)):
			totalmass += myData[j][6]
			if j < len(myData) - 1:
				for k in range(j+1,len(myData)):
				#print j,k
					Rmax = findRmax(myData[j][2], myData[j][3], myData[j][4], myData[k][2], myData[k][3], myData[k][4],Rmax)
		for i in range(len(myData)):
			#print i 
			r = findr(centre2[0],centre2[1],centre2[2],myData[i][2],myData[i][3],myData[i][4])
			#print Rmax,r,totalmass,myData[i][6],Velocity
			T = Tfric(Rmax,r,totalmass,myData[i][6], Velocity)
			print T
			string = "  %s" %(T)
			fil.write(string)
		fil.write("\n")
		f.close()
		
if __name__ == '__main__':
	fil = open("Times.txt", "w+")
	f = open("groups_all_boxes_at_z0p025_ascii_mergertimes.dat", "r")
	data = f.readlines()

	string = "SELECT GroupNumber, SubGroupNumber, CentreOfPotential_x, CentreOfPotential_y, CentreOfPotential_z, HalfMassRad_DM, Mass, Velocity_x, Velocity_y, Velocity_z FROM RefL0100N1504_Subhalo WHERE "
	realstring = []
	projection = []
        for line in data:
 		words = line.split()
                if words[0] == "#column":
	        	continue
                for i in range(5,len(words)):
			projection.append(words[2])
                        if words[i] != "-1":
	                        string += "(Mass = (SELECT MAX(Mass) FROM RefL0100N1504_Subhalo WHERE DescendantID = %s) AND DescendantID = %s) OR "%(words[i],words[i])
			if words[i] == "-1":
                                string = string[:-3]
				break
                                #print string
                realstring.append(string)
		string = "SELECT GroupNumber, SubGroupNumber, CentreOfPotential_x, CentreOfPotential_y, CentreOfPotential_z, HalfMassRad_DM, Mass, CentreOfMass_x, CentreOfMass_y, CentreOfMass_z, Velocity_x, Velocity_y, Velocity_z, GroupID FROM RefL0100N1504_Subhalo WHERE "
        con = sql.connect("kbx631", password = "YBGtrb40")
	indices = [41, '27', '28', 42, '27', '28', 70, '27', '28', 71, '27', '28', 111, '26', '28', 116, '26', '27', 118, '26', '28', 128, '26', '28', 129, '26', '26', 131, '26', '28', 136, '26', '28', 150, '26', '28', 163, '26', '28', 178, '26', '28', 230, '25', '27', 231, '25', '28', 237, '25', '27', 239, '25', '28', 241, '25', '27', 244, '25', '28', 269, '25', '27', 287, '25', '28', 288, '25', '28', 294, '25', '28', 298, '25', '28', 305, '25', '28', 320, '25', '28', 325, '25', '27', 328, '25', '28', 330, '25', '28', 356, '25', '26', 364, '24', '28', 367, '24', '26', 371, '24', '25', 377, '24', '27', 387, '24', '28', 389, '24', '28', 398, '24', '28', 399, '24', '27', 401, '24', '26', 402, '24', '28', 405, '24', '25', 406, '24', '27', 424, '24', '25', 429, '24', '26', 430, '24', '27', 431, '24', '28', 440, '24', '27', 442, '24', '25', 446, '24', '27', 449, '24', '28', 461, '24', '28', 464, '24', '26', 472, '24', '28', 479, '24', '28', 484, '24', '27', 485, '24', '27', 486, '24', '28', 508, '23', '28', 509, '23', '24', 511, '23', '26', 512, '23', '26', 521, '23', '28', 524, '23', '28', 529, '23', '28', 530, '23', '25', 532, '23', '26', 545, '23', '26', 546, '23', '26', 548, '23', '26', 551, '23', '27', 554, '23', '26', 555, '23', '25', 556, '23', '26', 561, '23', '26', 562, '23', '27', 563, '23', '27', 567, '23', '27', 571, '23', '25', 573, '23', '28', 574, '23', '28', 577, '23', '25', 581, '23', '25', 582, '23', '25', 587, '23', '25', 605, '23', '28', 607, '23', '26', 614, '23', '26', 619, '23', '26', 623, '23', '27', 624, '23', '27', 631, '23', '27', 633, '23', '26', 638, '23', '27', 640, '23', '25', 641, '23', '27', 643, '23', '27', 654, '23', '25', 681, '22', '28', 686, '22', '25', 688, '22', '28', 691, '22', '25', 695, '22', '27', 696, '22', '26', 701, '22', '26', 702, '22', '26', 706, '22', '25', 709, '22', '24', 716, '22', '25', 720, '22', '25', 722, '22', '25', 724, '22', '25', 729, '22', '25', 731, '22', '25', 732, '22', '25', 733, '22', '27', 737, '22', '25', 739, '22', '28', 740, '22', '26', 743, '22', '28', 744, '22', '27', 745, '22', '26', 746, '22', '28', 747, '22', '23', 751, '22', '25', 753, '22', '26', 754, '22', '26', 755, '22', '23', 757, '22', '28', 760, '22', '25', 763, '22', '25', 764, '22', '28', 766, '22', '27', 768, '22', '27', 770, '22', '27', 775, '22', '27', 785, '22', '27', 786, '22', '26', 793, '22', '28', 814, '22', '27', 816, '22', '28', 822, '22', '24', 826, '22', '27', 836, '22', '28', 838, '22', '26', 888, '21', '25', 889, '21', '26', 890, '21', '25', 893, '21', '24', 894, '21', '25', 898, '21', '23', 899, '21', '23', 900, '21', '27', 914, '21', '27', 916, '21', '23', 917, '21', '23', 918, '21', '28', 922, '21', '25', 924, '21', '24', 928, '21', '25', 934, '21', '25', 940, '21', '23', 941, '21', '25', 942, '21', '27', 952, '21', '23', 954, '21', '26', 956, '21', '26', 958, '21', '24', 960, '21', '26', 964, '21', '26', 965, '21', '23', 971, '21', '24', 973, '21', '27', 975, '21', '27', 979, '21', '26', 981, '21', '25', 982, '21', '24', 986, '21', '24', 990, '21', '26', 997, '21', '24', 1000, '21', '27', 1011, '21', '23', 1014, '21', '28', 1018, '21', '24', 1053, '20', '27', 1058, '20', '24', 1059, '20', '27', 1065, '20', '26', 1071, '20', '27', 1077, '20', '23', 1078, '20', '22', 1079, '20', '24', 1081, '20', '23', 1082, '20', '25', 1083, '20', '28', 1089, '20', '24', 1090, '20', '27', 1092, '20', '27', 1093, '20', '23', 1101, '20', '26', 1102, '20', '21', 1107, '20', '26', 1108, '20', '21', 1111, '20', '23', 1112, '20', '26', 1117, '20', '27', 1118, '20', '22', 1120, '20', '21', 1127, '20', '25', 1129, '20', '26', 1135, '20', '23', 1136, '20', '24', 1150, '20', '24', 1151, '20', '23', 1182, '19', '22', 1184, '19', '28', 1189, '19', '24', 1193, '19', '21', 1194, '19', '24', 1196, '19', '22', 1197, '19', '24', 1198, '19', '21', 1200, '19', '27', 1204, '19', '25', 1205, '19', '24', 1206, '19', '26', 1208, '19', '24', 1210, '19', '25', 1211, '19', '21', 1212, '19', '21', 1216, '19', '21', 1219, '19', '21', 1221, '19', '24', 1222, '19', '20', 1224, '19', '23', 1225, '19', '26', 1227, '19', '21', 1229, '19', '21', 1230, '19', '26', 1231, '19', '26', 1232, '19', '21', 1233, '19', '24', 1235, '19', '23', 1240, '19', '23', 1241, '19', '23', 1244, '19', '26', 1245, '19', '25', 1248, '19', '24', 1249, '19', '21', 1253, '19', '23', 1254, '19', '24', 1261, '19', '26', 1263, '19', '24', 1265, '19', '24', 1273, '19', '24']
        for i in range(0,len(indices),3):
	        myQuery = realstring[indices[i]]
                myData = sql.execute_query(con, myQuery)
                print myQuery
		#for j in range(len(myData)):
               		#print myData[j][13]
		string2 = "SELECT GroupMass, GroupCentreOfPotential_x, GroupCentreOfPotential_y, GroupCentreOfPotential_z, Group_R_Crit2500 FROM RefL0100N1504_FOF WHERE GroupID = %i" %myData[0][13]
		myData2 =  sql.execute_query(con, string2)
		myData2 = np.array([myData2])
		#print myData2[0][0]
		centre = np.array([myData2[0][1],myData2[0][2],myData2[0][3]])
		Findmass(myData[0][0], myData[0][1],  centre,i)
	fil.close()
