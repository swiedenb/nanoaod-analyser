import os

with open('run_2018_background_grid.txt', 'r') as fp:
    for line in fp:
        with open(line.split(':')[1].split('/')[-2] + '.json', 'w') as op:
            op.write('{\n')
            op.write(line.replace(',', ''))
            op.write('}')
        
