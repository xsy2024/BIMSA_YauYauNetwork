import yauyauSDE 
from Yau-Yau network import draw_graph
import os
import glob
import numpy as np
from scipy import stats

def linear_interplotion(t, time_step):  
      ts_new = []
      j = 0
      i = 0
      td = pd.DataFrame()
      while j <= 36:
          #print(j)
          if j in time_step:
              td[j] = t.iloc[:,i]
              i = i+1           
              j = j+0.25
          else: 
              td[j] = (t.iloc[:,i]*(j-time_step[i-1])+t.iloc[:,i-1]*(time_step[i]-j))/(time_step[i]-time_step[i-1])
              j = j+0.25
      return(td)

def main():
    folder_path = "./data"
    csv_files = glob.glob(os.path.join(folder_path, '*.csv'))
    data = {}
    data["time_step"] = np.array([0,0.5,1,1.5,2,4,6,8,10,12,16,20,24,30,36])
    for file in csv_files:
        print(f"正在读取文件: {file}")
        try:
            t = pd.read_csv(file, encoding='latin1')
        except UnicodeDecodeError:
            try:
                t = pd.read_csv(file, encoding='iso-8859-1')
            except UnicodeDecodeError:
                t = pd.read_csv(file, encoding='cp1252')
        t = t.set_index(t.columns[0],drop = True) 
        #t = np.log(t)
        t.loc["mean"] = [t.iloc[:,i].mean() for i in range(15)]
        t.loc["std"] = [t.iloc[:,i].std() for i in range(15)]
        t = t.apply(lambda x: np.where((x[:100] >=x.loc["mean"]+2*x.loc["std"]) & (x[:100]<= x.loc["mean"]-2*x.loc["std"]),x.loc["mean"],x[:100] ), axis = 0)
        _, opt = stats.boxcox(t.iloc[:100,:].values.ravel())
        t = (t**opt - 1)/opt
        data[f"{file}"] = t
    
    lambda_ = []
    for i in csv_files:
        #i = csv_files[0]
        t = data[i]
        _, opt = stats.boxcox(t.iloc[:100,:].values.ravel())
        print(opt)
        print(i)
        for j in range(1,15):
            dt = (t.iloc[:100,j]**opt - 1)/opt
            plt.plot(dt)
    
    # linear interplotion
    L_data = dict()
    for nm in csv_files[:3]:
        temp = linear_interplotion(data[nm], data["time_step"])
        L_data[nm] = temp
        L_data[f"d_{nm}"] = (temp.shift(-1, axis = "columns") - temp)/0.25
        L_data[f"d_{nm}"].iloc[:,-1] = L_data[f"d_{nm}"].iloc[:,-2]
      
    yysdes = yauyauSDE(csv_files[:3],L_data)
    yysdes.sub()
    yysdes.h()
    yysdes.L_conti(63, step = 0.25)
    yysdes.H_self(63)
    
    estimate = np.loadtxt("64_noise0.5.csv",delimiter=',')
    yysdes.L_conti(63, step = 0.001)
    #yysdes.w_deter("92",91)
    yysdes.w_new(estimate,63,"64")
  
    yysdes.Euler_gene(4, 0.001, 63)
    
    #X = np.array([yysdes.Total_data[j].iloc[91,:17] for j in yysdes.Total_name]).T
    yysdes.w_individual(state=estimate, name = "64_noise",no = 63)
    yysdes.w_deter( "64", no =63)
    #state = np.loadtxt("state_noise.csv",delimiter=',')
    X_0 = np.array([yysdes.Total_data[j].mean(0)[:4] for j in yysdes.Total_name]).T
    X_0 = yysdes.gene_.T
    yysdes.w(state=X_0, name="mean")

    G = nx.DiGraph()
    dot_num = len(name_list)
    time_step = list(range(400))
    ani = FuncAnimation(plt.gcf(), draw_graph, frames=time_step, interval = 50,fargs=(G,yysdes.Total_weight["64"],f'microbial 64 interactions',pos_param,name_list))
    ani.save("microbial_sample_ESP_3_63_new_yy.mp4",fps=20,dpi=500)

if __name__ == "__main__":
  
        main()
  
    
  
