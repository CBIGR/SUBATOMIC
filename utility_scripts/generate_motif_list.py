import argparse

def arg_parse():
    """
    Parses the arguments given in the command line
    """
    parser = argparse.ArgumentParser(description='Generates a list of preferable non-redundant motifs for a set of directed and a set of undirected edges.')
    parser.add_argument('-ui', help='Letters of the undirected networks seperated by _ . Example: P_G_H.')
    parser.add_argument('-di', help='Letters of the directed networks seperated by _ . Example: D_R_M.')
    parser.add_argument('-out', default=None, help='Path + name of the outfile. If not provided, result will be printed to console.')

    return parser.parse_args()


if __name__ == '__main__':
    """
    Generates a list of preferable non-redundant motifs for a set of directed and a set of undirected edges
    First argument: 
    """    
    
    # load arugments
    args = arg_parse()
    
    # process arguments
    try:
        ui = args.ui.split('_')
        di = args.di.split('_')

    except:
        print("please provide correct arguments for ui and di")
        print('-ui Letters of the undirected networks seperated by _ . Example: P_G_H.')
        print('-di Letters of the directed networks seperated by _ . Example: D_R_M.')
        exit()
   

    # calculate non-redundant motif selection based on input network letters

    # COM motifs
    COM = {}
    for i in ui:
        for j in ui:
            for k in ui:
                if i+j+k not in COM and i+k+j not in COM and j+i+k not in COM and k+i+j not in COM and j+k+i not in COM and k+j+i not in COM:
                    COM[i+j+k] = ''

    #print("COM", COM)

    # FFL motifs. DDD=Ddd=ddD. DDD motif alone is sufficient to represent all of the possibilities. However, all combinations of letters have to be kept
    FFL = {}
    for i in di:
        for j in di:
            for k in di:
                FFL[i+j+k] = ''
    #print("FFL", FFL)      
    

    # CIR motifs are either dDd or DdD. because it is circular, only one letter comination of all three letters is necessary
    CIR = {}
    for i in di:
        for j in di:
            for k in di:
                if i.upper()+j.lower()+k.upper() not in CIR and i.upper()+k.lower()+j.upper() not in CIR and j.upper()+i.lower()+k.upper() not in CIR and k.upper()+i.lower()+j.upper() not in CIR and j.upper()+k.lower()+i.upper() not in CIR and k.upper()+j.lower()+i.upper() not in CIR:
                    CIR[i.upper()+j.lower()+k.upper()] = ''
    
    #print("CIR", CIR)  

    

    # COP motifs UDD=DUd=ddU. If the position of the directed edges is switched, both motifs are redundant
    COP = {}
    for i in ui:
        for j in di:
            for k in di:
                if i+j+k not in COP and i+k+j not in COP:
                    COP[i+j+k] = ''

    #print("COP", COP)     
         
    # COR motifs, Udd=dUD=DDU. If the position of the directed edges is switched, both motifs are redundant
    COR = {}
    for i in di:
        for j in di:
            for k in ui:
                if i+j+k not in COR and j+i+k not in COR:
                    COR[i+j+k] = ''    
    #print("COR", COR)   
    
    # FBU motifs, UDd=UdD=DUD=dUd=DdU=dDU. We choose DUD. Changing the directed edges makes a difference, so no further motifs are redundant !
    FBU = {}
    for i in di:
        for j in ui:
            for k in di:
                FBU[i+j+k] = ''    
    #print("FBU", FBU)

     # FB2U motifs, UUD=UUd=UDU=UdU=DUU=dUU, We choose UUD. Changing the directed edges makes a difference, so no further motifs are redundant !
    FB2U = {}
    for i in ui:
        for j in ui:
            for k in di:
                FB2U[i+j+k] = ''    
    #print("FB2U", FB2U)

    # print output
    if args.out != None:
        outfile = open(args.out, 'w') 
    
    ALL = {**FFL, **COM, **CIR, **COP, **COR, **FBU, **FB2U} 

    for i in ALL:
        if args.out != None:
            outfile.write(i+'\n')
        else: 
            print(i)

    infile = open()

    
