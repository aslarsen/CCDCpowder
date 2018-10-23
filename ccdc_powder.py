import sys
from ccdc.descriptors import PowderPattern
from ccdc.io import CrystalReader
from os import remove
import numpy as np    
 
class PowderTrajectory:
    def __init__(self, pdb_trajectory):
        self._pdb_trajectory = pdb_trajectory
        self._author = "anders.sl@gmail.com"

    def run(self):
        self.read_pdbtrajectory(self._pdb_trajectory)

    def write_pdbfile(self, filename, pdbstring):
        f = open(filename, 'w')
        for line in pdbstring:
            f.write(line)
        f.close()

    def calc_powder_from_pdb(self, pdbfile):
        crystal = CrystalReader(pdbfile)[0]
        print crystal.spacegroup_symbol
        print 'angles', crystal.cell_angles[0], crystal.cell_angles[1], crystal.cell_angles[2]
        print 'length', crystal.cell_lengths[0], crystal.cell_lengths[1], crystal.cell_lengths[2]
        print 'volume', crystal.cell_volume

        pattern = PowderPattern.from_crystal(crystal)
        pattern.write_xye_file(('./'+pdbfile.replace('.pdb','.xye')))

    def calc_average_pattern(self, list_xye):
        patterns = []

        for xye in list_xye:
            data = []
            f = open(xye, 'r')
            for line in f:
                data.append(line.split())
            patterns.append(data)

        wavelength = patterns[0][0][0]
        thetas = []
        intensities = []
        uncertainties = []        
        variances = []

        for i in range(1,len(patterns[0])):
            theta = 0
            average_intensity = []
            average_uncertainty = []
            for j in range(0,len(patterns)):
                #print patterns[j][i][0], patterns[j][i][1], patterns[j][i][2]                            
                theta                = patterns[j][i][0]
                average_intensity.append(   float(patterns[j][i][1]))
                average_uncertainty.append( float(patterns[j][i][2]))
            variance = np.var(average_intensity)
            #print average_intensity, variance     
            average_intensity = np.mean(average_intensity)
            average_uncertainty = np.mean(average_uncertainty)
            #print theta, average_intensity, average_uncertainty, variance
            thetas.append(theta)
            intensities.append(average_intensity)
            uncertainties.append(average_uncertainty)
            variances.append(variance)

        xye_name = self._pdb_trajectory.replace('.pdb','_trajectory.xye')
        f = open(xye_name,'w')
        f.write(wavelength+'\n')
        for theta, intensity, uncertainty in zip(thetas, intensities, uncertainties):
            line = '{:6s} {:10.4f}  {:10.4f}'.format(theta, intensity, uncertainty) 
            #print line
            f.write(line+'\n')
        f.close()

        variance_name = self._pdb_trajectory.replace('.pdb','_variances.xye')
        f = open(variance_name,'w')
        for theta, variance in zip(thetas, variances):
            line = '{:6s} {:10.4f}'.format(theta, variance)
            f.write(line+'\n')
        f.close()

    def read_pdbtrajectory(self, pdb_trajectory):
        print 'Script by '+self._author

        with open(pdb_trajectory, 'r') as f:
            xye_filenames = []
            pdb_model = []
            model_number = 'none'
            first = True
            for line in f:                 
                if 'MODEL' in line:
                    model_number = line.split()[1]
                
                if 'ENDMDL' in line:
                    pdb_model.append(line)
                    print 'MODEL', model_number                    
                    temp_file = 'TEMP_'+model_number+'.pdb'
                    self.write_pdbfile( temp_file, pdb_model)     
                    self.calc_powder_from_pdb(temp_file)
                    xye_filenames.append(temp_file.replace('.pdb','.xye'))

                    remove( temp_file )

                    pdb_model = []
                else:
                    pdb_model.append(line)
            # merge xye files        
            self.calc_average_pattern(xye_filenames)

            print 'DONE!'
            

powder_trajectory = PowderTrajectory(sys.argv[1]) 
powder_trajectory.run()

