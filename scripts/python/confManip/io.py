'''
MODULE
  io
PURPOSE
  Analog of io_mod module of scf. Defines a class IO which
  provides generic interfaces for reading in scalars, vectors, 
  and matrices for several data types.  
COMMENT
  Because variables in python are not typed, only one input
  and one output function is required for io of scalar variables, 
  one for vectors, and one for matrices. In both input and output
  functions, the variable type is passed passed as a string 
  argument 'type', which can have values 'int', 'real', 'char' 
  and (for scalars) 'logical'.
'''

class IO:

    def __init__(self):
	self.comment = None

    def input_comment(self,file,comment=None):
        if not self.comment:
            self.comment = file.readline().strip()
        if comment:
            if comment == self.comment:
                self.comment = None
                return 1
            else:
                return 0
        else:
           self.comment = None
           return 1
            
    # Scalar Variables
        
    def input_var(self,file,type,comment=None,f='A'):
        '''
        PURPOSE
          Read and return a scalar variable of specified type 
          from a specified file 
        ARGUMENTS
          file - file object (must be opened for reading)
          type - string, 'int', 'real', 'char', or 'logic'
          f    - 'A' -> comment string on line above data
               - 'N' -> no comment string
        '''
        if f == 'A':
        #    self.comment = file.readline().strip()
	    if not self.input_comment(file,comment):
	        return None
        data = file.readline().strip()
        if type == 'int':
            return int(data)
        elif type == 'real':
            return float(data)
        elif type == 'char':
            return strip_quotes(data)
        elif type == 'logic':
            if data in ['T','true', 'TRUE', 'True']:
                return 1
            elif data in ['F','false','FALSE','FALSE']:
                return 0
            else:
                raise "Invalid logical variable:", data
        else:
            raise 'Illegal type in input_var'
    
    # Vectors
        
    def input_vec(self,file,type,n=None, comment=None, s='R', f='A'):
        '''
        PURPOSE
          Read and return a vector of n variables of specified 
          type from file 
        ARGUMENTS
          file - file object (must be opened for reading)
          type - string, = 'int', 'real', 'char', or 'logic'
          n    - number of elements in vector
          f    - 'A' -> comment on line above data
               - 'N' -> no comment 
          s    - 'R' -> row vector (one line)
               - 'C' -> column vector (n lines)
        '''
        if f == 'A':
            #self.comment = file.readline().strip()
	    if not self.input_comment(file,comment):
	        return None
        if s == 'R':   # Row Vector
            data = file.readline().split()
            if n:
                data = data[:n]
        elif s == 'C': # Column Vector
            if not n: raise 'No n value in input_int_vec for column vector'
            data = []
            for i in range(n):
                data.append( file.readline().split()[0] )
        if type == 'int':
            return [ int(x) for x in data ]
        elif type == 'real':
            return [ float(x) for x in data ]
        elif type == 'char':
            return [ strip_quotes(x) for x in data ]
        else:
            raise 'Illegal type in input_vec'
    
    # Matrices
        
    def input_mat(self,file,type,m,n=None,comment=None,s=None,f='A'):
        '''
        PURPOSE
          Read and return a m x n matrix of variables of specified 
          type from file 
        ARGUMENTS
          file - file object (must be opened for reading)
          type - string, = 'int' or 'real'
          m    - number of rows in matrix
          n    - number of columns (set to m by default)
          f    - 'A' -> comment on line above data
               - 'N' -> no comment 
          s    - 'N' -> full matrix
               - 'S' -> symmetric
               - 'L' -> symmetric lower triangular (0 diagonals)
        '''
        if f == 'A':
	    if not self.input_comment(file,comment):
	        return None
        # Emulate initialization of m x n array
        if not n:
           n = m
        data = []
        for i in range(m):
            data.append([])
            for j in range(n):
                if type == 'int':
                    data[i].append(0) 
                elif type == 'real':
                    data[i].append(0.0) 
        # Read matrix
        if s == 'N' or s == 'S':
            min = 0
        elif s == 'L':
            min = 1
        else:
            raise 'Invalid style in input_mat: s=' + str(s)
        for i in range(min,m) :
            line = file.readline().split()
            if s == 'N':
                line = line[:n]
            elif s == 'S':
                line = line[:i+1]
            elif s == 'L':
                line = line[:i+1]
            for j in range(len(line)):
                if type == 'int':
                    datum = int(line[j])
                elif type == 'real':
                    datum = float(line[j])
                else:
                    raise 'Illegal type in input_mat'
                data[i][j] = datum
                if s == 'S' or s == 'L':
                    data[j][i] = datum
        return data
    
    
    def format_var(self,type,data):
        '''
        PURPOSE
            Returns a formatted string representation of data,
    	formatted appropriately for specified data type
        '''
        if type == 'int':
            return '%15d' % data
        elif type == 'real':
            return '%15.7E' % data
        elif type == 'char':
    	    data = "'" + data + "'" 
	    return "%15s" % data 
        elif type == 'logic':
            if data:
                data = 'T'
            else:
                data = 'F'
            return "%15s" % data 
        else:
            raise 'Illegal type in format_var'
      
    def output_comment(self,file, comment, n=20, f='w'):
        comment = comment.strip()
        comment = comment.ljust(n)
        if f == 'w':
            comment = file.write(comment + "\n")
        else:
            comment = file.write(comment + "    ")
    
    def output_var(self, file, type, data, comment, f='B'):
        '''
        PURPOSE
          Output a scalar variable of specified type to a file
        ARGUMENTS
          file - file object (must be opened for reading)
          type - string, 'int', 'real', 'char', or 'logic'
          f    - 'A' -> comment string on line above data
               - 'B' -> comment string on the same line as data
               - 'N' -> no comment string
        '''
        if f == 'A':
            self.output_comment(file,comment)
        elif f == 'B':
            self.output_comment(file,comment,20,'n')
        file.write(self.format_var(type,data) + "\n")
    
    def output_vec(self,file, type, data, n=None, comment=None, s='R', f='A'):
        '''
        PURPOSE
          Output a vector of n variables of specified type to a file 
        ARGUMENTS
          file - file object (must be opened for reading)
          type - string, = 'int', 'real', 'char', or 'logic'
          n    - number of elements in vector
          f    - 'A' -> comment on line above data
               - 'N' -> no comment 
          s    - 'R' -> row vector (one line)
               - 'C' -> column vector (n lines)
        '''
	if not n:
  	    n = len(data)
        if f == 'A':
            self.output_comment(file,comment)
        if f == 'B':
            self.output_comment(file,comment,20,'n')
        for i in range(n):
            file.write(self.format_var(type,data[i]))
            if s == 'C':
                file.write("\n")
        if s == 'R':
            file.write("\n")
    
    
    # Matrices
        
    def output_mat(self,file,type,data,m=None,n=None,comment=None,s='L',f='A'):
        '''
        PURPOSE
          Output a m x n matrix of variables of specified type to file
        ARGUMENTS
          file - file object (must be opened for reading)
          type - string, = 'int' or 'real'
          m    - number of rows in matrix
          n    - number of columns (set to m by default)
          f    - 'A' -> comment on line above data
               - 'N' -> no comment 
          s    - 'N' -> full matrix
               - 'S' -> symmetric
               - 'L' -> symmetric lower triangular (0 diagonals)
        '''
	if not m:
           m = len(data)
        if not n:
           n = m
        if f == 'A':
            self.output_comment(file,comment)
        # Read matrix
        if s == 'N' or s == 'S':
            min = 0
        elif s == 'L':
            min =1
        for i in range(min,m) :
            if s == 'N':
                max = n
            elif s == 'S':
                max = i+1
            elif s == 'L':
                max = i
            vec = data[i][:max]
            self.output_vec(file, type, vec, max, None, s='R', f='N')


def strip_quotes(q_string):
    '''
    Strip initial and final quotes marks from string q_string,
    if present. Return stripped string. 
    '''
    q_string = q_string.strip()
    if q_string[0] == "'" and q_string[-1] == "'":
        return q_string[1:-1]
    elif q_string[0] == '"' and q_string[-1] == '"':
        return q_string[1:-1]
    else:
        return q_string

