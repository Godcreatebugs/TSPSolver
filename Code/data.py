import random

with open("input300_10.txt","wb") as the_file:
    the_file.write('NAME : input300_10.txt \n')
    the_file.write('COMMENT : \n')
    the_file.write('TYPE : 2D Metric Graph\n')
    the_file.write('DIMENSION : 300 \n')
    the_file.write('EDGE_WEIGHT_TYPE : Metric 2D\n')
    the_file.write('NODE_COORD_SECTION\n')
    for i in range(300):
        the_file.write(str(i+1))
        the_file.write(" ")
        the_file.write(str(random.randint(1,200)))
        the_file.write(" ")
        the_file.write(str(random.randint(1,200)))
        the_file.write('\n')
the_file.close()
