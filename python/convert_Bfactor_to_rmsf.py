import math

def get_all_b_factors(in_filename, out_filename):

  input_handle = open(in_filename, 'r')

  output_handle = open(out_filename, 'w')
  
  output_handle.write('converted_rmsf\n')

  for line in input_handle.readlines():
    if len(line) == 81 and line[0:4] == 'ATOM':
      line = str((3*float(line[61:66])/(8*math.pi**2))**(1./2)) + '\n'
      output_handle.write(line)

  output_handle.close()

  return 0
  
def get_ca_b_factors(in_filename, out_filename):

  input_handle = open(in_filename, 'r')

  output_handle = open(out_filename, 'w')
  
  output_handle.write('converted_rmsf\n')

  for line in input_handle.readlines():
    if len(line) == 81 and line[0:4] == 'ATOM' and line[13:15] == 'CA':
      line = str((3*float(line[61:66])/(8*math.pi**2))**(1./2)) + '\n'
      output_handle.write(line)

  output_handle.close()

  return 0



