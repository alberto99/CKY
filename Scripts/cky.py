
class SetupCKY(object):
    def __init__(self,levels=2):
        self.levels = levels
        self.trajectories = []
        self.trajectory_index = {}
        self.start_index()
        self.combination_rules()
        self.possible_simulation_fragments()
	self.cycles=1000
        #self.CKY_combinations()

    def start_index(self):
        self.index = {}
        self.index[0] = []
        self.index[0].append(0)

    def combination_rules(self):
        #Look for possible combinations of CKY elements
        self.comb = {}
        self.combinations = []
        for lev in range(1,self.levels+1):
            for x in range(1,lev+1):
                for y in range(0,x+1):
                    self.process(lev,x,y)
        #add the last one where you just do a simulation no CKY involved
        #self.process(self.levels,self.levels,0)

    def possible_simulation_fragments(self):
        #Transform all possible combinations into which elements are being combined
        #remove symmetry, but keep the different paths
        #eg 12+3 == 3+12; but 1+23 != 12+3

        #all units to work with in parse tree
        states = []
        for i in range(1,self.levels+1):
            for j in range(1,self.levels-i+2):
                state = j
                for k in range(1,i):
                    state = '%s%s' % (state,j+k)
                states.append(state)
        self.states = states

    def CKY_combinations(self):
        #Combining the possibilities
        #Into the CKY trajectories we are going to carry out
        for state in self.states:
            a = len("%s" % state)
            options = self.comb[a]
            self.combine(state,a,options)

    def indexing(self,state,el):
        state = '%s' % state
        try:
            self.index[state].append(el)
        except:
            self.index[state] = []
            self.index[state].append(el)

    def combine(self,state,a,options):
        s = "%s" % state
        #easiest case, no combination
        if a == 1:
            el = CKYelement(s, state, 0)
            el.txt = 'Single sse trajectory: %s' % s
            el.parent1 = None
            el.index = len(self.trajectories)
            self.trajectories.append(el)
            p1 = '%s_%s' % (None,None)
            p2 = '%s_%s' % (None,None)
            self.trajectory_index[(state,0,p1,p2)] = self.trajectories.index(el)
            self.indexing(state,el)
            return
        #For two elements there are no possible combinations except 1 + 1
        if a == 2:
            el = CKYelement(state, s[0],s[1])
            el.txt = 'State %s produced from %s and %s' % (s,s[0],s[1])
            el.parent_index1 = self.trajectory_index[(int(s[0]),0,'None_None','None_None')]
            el.parent_index2 = self.trajectory_index[(int(s[1]),0,'None_None','None_None')]
            el.index = len(self.trajectories)
            self.trajectories.append(el)
            p1 = '%s_%s' % (s[0],None)
            p2 = '%s_%s' % (s[1],None)
            self.trajectory_index[(int(s[0]),int(s[1]),p1,p2)] = self.trajectories.index(el)
            self.indexing(state,el)
            return
        #More elements. Now certain elements can have different origins. Make a trajectory from each
        #possibility. E.g. 123 can come from 12+3 and from 1+23
        for (x,y) in options:
            e1 = s[0:x]
            e2 = s[x:]
            if not e2:
                e2 = 0
            for t_e1 in self.index[e1]:
                for t_e2 in self.index[e2]:
                    self.combine_higher(s,t_e1,t_e2,e1,e2)

            if (x == y) or (y==0):
                continue
            #otherwise, there is another possibility:
            e1 = s[0:y]
            e2 = s[y:]
            for t_e1 in self.index[e1]:
                for t_e2 in self.index[e2]:
                    self.combine_higher(s,t_e1,t_e2,e1,e2)

    def combine_higher(self,s,t_e1,t_e2,e1,e2):
                    try:
                        el = CKYelement(s,t_e1.state,t_e2.state)
                        el.txt = 'State %s produced from %s and %s,coming from (%s _ %s) and (%s _ %s)' % (s,t_e1.state,t_e2.state,t_e1.parent1,t_e1.parent2,t_e2.parent1,t_e2.parent2)
                        el.name = 'sse%s' % s
                        el.parent_index1 = self.trajectories.index(t_e1)
                        el.parent_index2 = self.trajectories.index(t_e2)
                        el.index = len(self.trajectories)
                        self.trajectories.append(el)
                        p1 = '%s_%s' % (t_e1.parent1,t_e1.parent2)
                        p2 = '%s_%s' % (t_e2.parent1,t_e2.parent2)
                        self.trajectory_index[(int(e1),int(e2),p1,p2)] = self.trajectories.index(el)
                        self.indexing(s,el)
                    except:
                        #Avoid endless loop when we do the non CKY trajectory

                        #Additionally, don't add if we already created
                        try:
                            index = self.trajectory_index[(int(e1),int(e2),'None_None','None_None')]
                        except:
                            el = CKYelement(s,e1,e2)
                            el.txt = 'State %s produced from %s and %s' % (s,e1,e2)
                            el.name = 'sse%s' % s
                            el.index = len(self.trajectories)
                            self.trajectories.append(el)
                            p1 = '%s_%s' % (None,None)
                            p2 = '%s_%s' % (None,None)
                            self.trajectory_index[(int(e1),int(e2),p1,p2)] = self.trajectories.index(el)
                            self.indexing(s,el)

    def process(self,lev,x,y):
        add = x + y
        if add != lev:
           return 
        if 1 < x < self.levels:
            if y == 0:
                return
        self.combinations.append( (lev,x,y) )
        try:
            self.comb[lev].append( (x,y) )
        except:
            self.comb[lev] = []
            self.comb[lev].append( (x,y) )


class CKYelement(object):
    def __init__(self,state,sse1,sse2):
        self.sse1 = sse1
        self.sse2 = sse2
        #print "CKY", state
        self.state = state
        self.name = 'sse%s' % state
        self.parent1 = sse1
        self.parent2 = sse2
        if self.parent1 == 0:
            self.parent1 = None
        if self.parent2 == 0:
            self.parent2 = None
        self.parent_index1 = None
        self.parent_index2 = None



def main():
    CKYtrajs = SetupCKY(levels=4)
    print CKYtrajs.index.keys()
    for e in CKYtrajs.trajectories:
        #print e.txt
        print e.name,e.state
        try:
            #print e.parent1,e.parent2,e.parent_index1,e.parent_index2
            p1 = e.parent_index1
            p2 = e.parent_index2
            t1 = CKYtrajs.trajectories[p1]
            t2 = CKYtrajs.trajectories[p2]
            #print t1.parent1,t1.parent2
            #print t2.parent1,t2.parent2
        except:
            #print e.parent1,e.parent2
            pass

    print len(CKYtrajs.trajectories)



if __name__ == '__main__':
    main()
