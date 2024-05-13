from fractions import Fraction
from typing import Union

def max_min_index(array,is_max,defult_border = None):
    if len(array) <= 0:
        return None
    border = defult_border
    result = None
    for i in range(len(array)):
        if array[i] == None:
            continue
        if border == None or ((array[i] > border) if is_max else (array[i] < border)):
            result = i
            border = array[i]
    return result
    
def max_index(array:list,minnum:int = None) -> Union[int,None]:
    '''获取数组中大于 `minnum` 的最大值所在位置'''
    return max_min_index(array,True,minnum)

def min_index(array:list,maxnum:int = None) -> Union[int,None]:
    '''获取数组中小于 `maxnum` 的最小值所在位置'''
    return max_min_index(array,False,maxnum)

class Simplex():
    '''单纯形法复合求解器'''
    def __init__(self) -> None:
        self.A = []
        self.B = []
        self.C = []
        self.baseIndexes = []
        self.last_enbase = None
        self.last_unbase = None
        self.Cb = []
        self.names = []

    def index_minB_below_zero(self):
        '''返回 B 向量非正最小值位置'''
        return min_index(self.B,0)
    
    def check_number(self):
        '''返回检验数数组'''
        result = self.C.copy()
        for i in range(len(self.B)):
            for j in range(len(self.C)):
                result[j] -= self.A[i][j] * self.Cb[i]
        return result
    
    def get_theta(self,enbase_index):
        '''获取 θ 列'''
        theta = []
        for i in range(len(self.B)):
            if (aij := self.A[i][enbase_index]) > 0:
                theta.append(self.B[i]/aij)
            else:
                theta.append(None)
        return theta

    def solveStep(self):
        '''使用单纯形法解决单步问题'''
        if index := self.index_minB_below_zero():
            return self.__solveDualStep(minB_index=index)
        else:
            return self.__solveNormalStep()

    def get_current_solution(self):
        '''获取当前解'''
        sol = []
        for j in range(len(self.C)):
            if j in self.baseIndexes:
                sol.append(self.B[self.baseIndexes.index(j)])
            else:
                sol.append(0)
        return sol

    def get_current_solution_str(self):
        '''获取当前解字符串'''
        sol = self.get_current_solution()
        return f'[{str.join(', ',(str(v) for v in sol))}]'
    
    def __check_inf_solution(self,sigma) -> bool:
        '''无穷解检查'''
        for i in range(len(self.C)):
            if i in self.baseIndexes:
                continue
            if sigma[i] == 0:
                print(f" [!] {i} 列非基变量检验数为 0，该线性规划问题存在无穷解！")
                solution1 = self.get_current_solution_str()
                enbase_index = i 
                theta = self.get_theta(i)
                unbase_index = min_index(theta,None)
                self.print_table(sigma,theta,enbase_index,unbase_index)
                self.last_enbase = enbase_index
                self.last_unbase = unbase_index
                self.__solveNormalStep(True)
                solution2 = self.get_current_solution_str()
                print(f"该问题的两个解为：{solution1} 和 {solution2}")
                return True
        else:
            return False
            
    def __solveNormalStep(self,in_infsol_search = False):
        '''使用正常单纯形法解决该步骤'''
        # 初等行列变换
        unbase_index = self.last_unbase
        enbase_index = self.last_enbase
        if self.last_enbase != None and self.last_unbase != None:
            div = self.A[unbase_index][enbase_index]
            for j in range(len(self.C)):
                self.A[unbase_index][j] /= div
            self.B[unbase_index] /= div
            for i in range(len(self.B)):
                if i == unbase_index:
                    continue
                times = self.A[i][enbase_index]
                for j in range(len(self.C)):
                    self.A[i][j] -= times * self.A[unbase_index][j]
                self.B[i] -= times * self.B[unbase_index]
        # 寻找进基出基变量
        sigma = self.check_number()
        if (enbase_index := max_index(sigma,0)) == None:
            if not in_infsol_search and self.__check_inf_solution(sigma): # 无穷解检查
                return True
            self.print_table(sigma,None,enbase_index,unbase_index)
            if not in_infsol_search:
                print(f"该问题的解为：{self.get_current_solution_str()}")
            return True
        theta = self.get_theta(enbase_index)
        unbase_index = min_index(theta,None)
        self.print_table(sigma,theta,enbase_index,unbase_index)
        if unbase_index == None:
            print("在该表中找不到出基变量，该线性规划问题存在无界解！")
            return True
        self.baseIndexes[unbase_index] = enbase_index # 转移基变量
        self.Cb[unbase_index] = self.C[enbase_index]
        self.last_enbase = enbase_index
        self.last_unbase = unbase_index
        return False

    def __solveDualStep(self,minB_index):
        '''使用对偶单纯形法解决该步骤'''
        raise NotImplementedError('尚未实现')

    def fractionlize(self):
        '''分数实例化矩阵'''
        for i in range(len(self.B)):
            if not isinstance(self.B[i],Fraction):
                self.B[i] = Fraction(self.B[i])
        for j in range(len(self.C)):
            if not isinstance(self.C[j],Fraction):
                self.C[j] = Fraction(self.C[j])
        for i in range(len(self.B)):
            for j in range(len(self.C)):
                if not isinstance(self.A[i][j],Fraction):
                    self.A[i][j] = Fraction(self.A[i][j])
   
    def formalize(self) -> None:
        '''将 LP 问题转化为标准型与求解器初始化'''
        self.Cb = []
        self.names = []
        for i in range(len(self.C)):
            self.names.append(f'x{i+1}')
        for i in range(len(self.B)):
            self.Cb.append(0)
            self.names.append(f'x{i+len(self.C)+1}')
            for j in range(len(self.B)):
                if i == j:
                    self.baseIndexes.append(len(self.C) + i)
                    self.A[i].append(Fraction(1))
                else:
                    self.A[i].append(Fraction(0))
        for i in self.B:
            self.C.append(0)

    def slove(self) -> None:
        '''解决线性规划问题'''
        self.fractionlize()
        self.formalize()
        i = 1
        while(True):
            print(f'第 {i} 次迭代:')
            if self.__solveNormalStep():
                break
            i += 1
        maxz = 0
        sol = self.get_current_solution()
        for j in range(len(self.C)):
            maxz += self.C[j] * sol[j]
        print(f'最优值 Z* = {maxz}')
    
    def print_table(self,sigma,theta,enbase_index,unbase_index):
        '''打印单纯形表'''
        # row 1
        print('|=======< C >=======>\t',end = '')
        for c in self.C:
            print(c,end = '\t')
        print()
        # row 2
        print('Cb\tBase\tb\t',end ='')
        for i in range(len(self.C)):
            print(self.names[i],end = '\t')
        if theta:
            print('θ')
        else:
            print()
        # row b1...bn
        for i in range(len(self.B)):
            print(self.Cb[i],end = '\t')
            print(self.names[self.baseIndexes[i]],end = '\t')
            print(self.B[i],end = '\t')
            for j in range(len(self.C)):
                if i == unbase_index and j == enbase_index:
                    print(f"[{self.A[i][j]}]",end = '\t')
                else:
                    print(self.A[i][j],end = '\t')
            if theta:
                print(theta[i] if theta[i] else '-')
            else:
                print()
        # row sigma
        print('|=======< σ >=======>\t',end = '')
        for j in range(len(self.C)):
            print(sigma[j],end = '\t')
        print('\n')
    
    
#lp = Simplex()
#lp.C = [3,3]
#lp.B = [12,16,15]
#lp.A = [[2,2],[4,0],[0,5]]
#lp.slove()
