import numpy as np
import eagleSqlTools as sql
import matplotlib.pyplot as plt
from sphviewer.tools import QuickView
from read_eagle import EagleSnapshot
import h5py

parttype = 4

if parttype == 0:
	parttype_string = 'Gas'
	max_hsml=1e5   # any very large number
elif parttype == 1:
	parttype_string = 'DarkMatter'
	max_hsml=1e5
elif parttype == 4:
	parttype_string = 'Stars'
	max_hsml=0.001  # a very small number; no smoothing for stars

def  read_header():
        f = h5py.File('/disks/eagle/L0100N1504/REFERENCE/data/snapshot_027_z000p101/snap_027_z000p101.0.hdf5','r')
	a = f['Header'].attrs.get('Time')
	h = f['Header'].attrs.get('HubbleParam')
	boxsize = f['Header'].attrs.get('Boxsize')
	f.close()
	return a,h, boxsize

def defineprojection(projection):
	print projection
	if projection == 'posx':
		p = 90
		t = 0
	if projection == 'posy':
		p = 0
		t = 270
	if projection == 'posz':
		p = 0
		t = 0  
	if projection == 'negx':
		p = 270
		t = 0
	if projection == 'negy':
		p = 0
		t = 90
	if projection == 'negz':
		p = 180
		t = 0
	return p,t

def find_galaxy(projection,galx,galy,galz):
	if projection == 'posx':
		x1 = -galz
		y1 = galy
		x2 = -galx
		y2 = galy
		x3 = galy
		y3 = galz
		l1 = "-z [Mpc]"
		l2 = "y [Mpc]"
		l3 = "-x [Mpc]"
		l4 = "y [Mpc]"
		l5 = "y [Mpc]"
		l6 = "z [Mpc]"
	if projection == 'posy':
		x1 = galx
		y1 = -galz
		x2 = -galy
		y2 = -galz
		x3 = galx
		y3 = galy
                l1 = "x [Mpc]"
		l2 = "-z [Mpc]"
		l3 = "-y [Mpc]"
		l4 = "-z [Mpc]"
		l5 = "x [Mpc]"
		l6 = "y [Mpc]"

	if projection == 'posz':
	        x1 = galx
		y1 = galy
		x2 = -galz
		y2 = galy
		x3 = galx
		y3 = galz
		l1 = "x [Mpc]"
		l2 = "y [Mpc]"
		l3 = "-z [Mpc]"
		l4 = "y [Mpc]"
		l5 = "x [Mpc]"
		l6 = "z [Mpc]"

	if projection == 'negx':
		x1 = galz
		y1 = galy
		x2 = galx
		y2 = galy
		x3 = -galy
		y3 = galz
		l1 = "z [Mpc]"
		l2 = "y [Mpc]"
		l3 = "x [Mpc]"
		l4 = "y [Mpc]"
		l5 = "-y [Mpc]"
		l6 = "z [Mpc]"

	if projection == 'negy':
		x1 = galx
		y1 = galz
		x2 = galy
		y2 = galz
		x3 = galx
		y3 = -galy
                l1 = "x [Mpc]"
		l2 = "z [Mpc]"
		l3 = "y [Mpc]"
		l4 = "z [Mpc]"
		l5 = "x [Mpc]"
		l6 = "-y [Mpc]"

	if projection == 'negz':
		x1 = -galx
		y1 = galy
		x2 = galz
		y2 = galy
		x3 = -galx
		y3 = galz
                l1 = "-x [Mpc]"
		l2 = "y [Mpc]"
		l3 = "z [Mpc]"
		l4 = "y [Mpc]"
		l5 = "-x [Mpc]"
		l6 = "z [Mpc]"

	return x1,y1,x2,y2,x3,y3,l1,l2,l3,l4,l5,l6
	
class PhaseDiagram_ReadEagle:

        def __init__(self,gn, sgn, centre, load_region_length=2):
		self.a, self.h, self.boxsize = read_header()
		self.gas = self.read_galaxy(parttype, gn, sgn, centre, load_region_length)

	def read_galaxy(self, itype, gn, sgn, centre, load_region_length):

		p,t = defineprojection(projection[indices[n]])
		data = {}
		
		eagle_data = EagleSnapshot('/disks/eagle/L0100N1504/REFERENCE/data/snapshot_027_z000p101/snap_027_z000p101.0.hdf5')
		
		centre *= self.h
			
		region = np.array([(centre[0]-0.5*load_region_length), (centre[0]+0.5*load_region_length), (centre[1]-0.5*load_region_length), (centre[1]+0.5*load_region_length), (centre[2]-0.5*load_region_length),(centre[2]+0.5*load_region_length)])
		
		eagle_data.select_region(*region)
		f = h5py.File('/disks/eagle/L0100N1504/REFERENCE/data/snapshot_027_z000p101/snap_027_z000p101.0.hdf5', 'r')

		constants = f['Constants']
	        Mpc = constants.attrs['CM_PER_MPC']
	        SMass = constants.attrs["SOLAR_MASS"]
		for att in ['GroupNumber','Coordinates','Mass']:
			if parttype == 1 and att == 'Mass':
				data[att] = np.ones(len(data['GroupNumber'])) * SMass * 9.7e6
				continue
			tmp = eagle_data.read_dataset(itype, att)
			cgs = f['PartType%i/%s'%(itype,att)].attrs.get('CGSConversionFactor')
			aexp = f['PartType%i/%s'%(itype,att)].attrs.get('aexp-scale-exponent')
			hexp = f['PartType%i/%s'%(itype,att)].attrs.get('h-scale-exponent')
			data[att] = np.multiply(tmp, cgs * self.a**aexp * self.h**hexp, dtype = 'f8')
			
			if att == 'Coordinates':
				centre2 = np.multiply(centre, cgs / Mpc * self.a**aexp *self.h**hexp)
                                galx = np.zeros(len(myData))
			        galy = np.zeros(len(myData))
			        galz = np.zeros(len(myData))
			        j = 1
		                print myData[1][2], cgs, Mpc, self.a, aexp, self.h, hexp
				print np.multiply(myData[1][2],cgs / Mpc * self.a**aexp * self.h**hexp)
				for i in range(1,len(myData)):  #Only if first galaxy in list is centre
					galx[j] = np.multiply(myData[i][2],cgs / Mpc * self.a**aexp * self.h**hexp *self.h) - centre2[0]
					galy[j] = np.multiply(myData[i][3],cgs / Mpc * self.a**aexp * self.h**hexp *self.h) - centre2[1]
                                	galz[j] = np.multiply(myData[i][4],cgs / Mpc * self.a**aexp * self.h**hexp *self.h) - centre2[2]
                                	j += 1
	        		print galx,galy,galz
				x1,y1,x2,y2,x3,y3,l1,l2,l3,l4,l5,l6 = find_galaxy(projection[indices[n]],galx,galy,galz)

		f.close()

		mask = (data['GroupNumber'] == gn)
		for att in data.keys():
			data[att] = data[att][mask]	
		for i in range(3):
			if i == 0:
				stringtitle = parttype_string + " for p=0 and t=0" + " at snapshot 27 with redshift z = 0.1"
				x = x1
				y = y1
				xlabel = l1
				ylabel = l2
			if i == 1: 
				p = p + 90
				stringtitle = parttype_string + " for p=90 and t=0" + " at snapshot 27 with redshift z = 0.1"
				x = x2
				y = y2
				xlabel = l3
				ylabel = l4
			if i == 2:
				stringtitle = parttype_string + " for p=0 and t=90" + " at snapshot 27 with redshift z = 0.1"
				t = t + 90
				p = p - 90
				x = x3
				y = y3
				xlabel = l5
				ylabel = l6
			qv_parallel = QuickView(data['Coordinates']/Mpc, data['Mass']/SMass, r = 'infinity', plot = False,p=p, t=t,  x=centre2[0], y = centre2[1], z = centre2[2],max_hsml = max_hsml, extent=[-0.5,0.5,-0.5,0.5])

	
			plt.imshow(qv_parallel.get_image(), extent=qv_parallel.get_extent(), cmap ='viridis', origin='lower')
			plt.colorbar()
			color = ['r', 'k', 'm', 'w']
			plt.title(stringtitle)
			plt.xlabel(xlabel)
			print p,t,x,y
			plt.ylabel(ylabel)
			plt.xlim(-0.21,0.21)
			plt.ylim(-0.21,0.21)
			#plt.clim(5,)
			#plt.scatter(x, y, facecolors = 'none', color = 'r')
                        marker = ['o', 's', '^', 'D']
                        for j in range(len(x)):
	                        plt.scatter(x[j], y[j],marker = marker[j], facecolors = 'none', color = 'r', label = myData[j][5])
                        plt.legend(bbox_to_anchor=(-0.3,-0.15), loc=3, prop={'size': 8})

			plt.show()
		
if __name__ == '__main__':

	f = open("groups_all_boxes_at_z0p025_ascii.dat", "r")
	#f = open("groups_all_boxes_at_z0p025_ascii_mergertimes.dat", "r")
	data = f.readlines()

	string = "SELECT GroupNumber, SubGroupNumber, CentreOfPotential_x, CentreOfPotential_y ,CentreOfPotential_z, GalaxyID FROM RefL0100N1504_Subhalo WHERE "
	realstring = []
	projection = []
	for line in data:
		words = line.split()
		if words[0] == "#column":
			continue
		for i in range(3,len(words)):
		#for i in range(5,len(words)):
			projection.append(words[2])
			if words[i] != "-1":
				string += "GalaxyID = %s or "%(words[i])
			if words[i] == "-1":
				string = string[:-3]
				break
		realstring.append(string)
		string = "SELECT GroupNumber, SubGroupNumber, CentreOfPotential_x, CentreOfPotential_y ,CentreOfPotential_z, GalaxyID FROM RefL0100N1504_Subhalo WHERE "
	con = sql.connect("kbx631", password = "YBGtrb40")
	indices = [11,13,19,24,49,64,69]
	#indices = [42]
	n = 4
	for i in range(3,len(indices)):
		myQuery = realstring[indices[i]]
		print myQuery
		myData = sql.execute_query(con, myQuery)
		centre = np.array([myData[0][2],myData[0][3],myData[0][4]])
        	PhaseDiagram_ReadEagle(myData[0][0], myData[0][1],  centre)
		n += 1
#	xp = PhaseDiagram_ReadEagle(941030,0, centrep)
#	y = PhaseDiagram_ReadEagle(1846,2, centre2)
