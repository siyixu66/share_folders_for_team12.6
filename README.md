# share_folders_for_team12.6
天体测量及其应用小组作业

selcon.f里的236行换成JPLEPH的绝对路径      NAMFIL='/Users/jcyuan/Desktop/lunar/JPLEPH'
然后编译一下就可生成可执行文件了

test.f90为主程序，可以自行设置开始年份和结束年份
selcon.f里是调用星历及常量的子函数
compute.f是sofa程序库里关于时间转换的函数

目前计算的初亏时间：
2025 年           3 月          14 日           5 :          10 :   18.0119991   
2025 年           9 月           7 日          16 :          27 :   39.9900017    

和NASA给出的时、分一致，NASA没有给出秒：https://eclipse.gsfc.nasa.gov/JLEX/JLEX-AS.html
计算了2019年的月食，和紫金山天文台给出的初亏时间对比，时、分一致，秒不一样

# 12.12
添加了新的主程序main.f90,主要加入了顾及地月光行时的代码块，部分变量有改动

# 12.15
上传了更新后的lunar_eclipse2.f90,该程序可以计算完整的月食过程，并且计算食分，判断是月偏食还是月全食

# 1.2
上传了lunar_eclipse3.f90,该程序增加了地月光行时迭代
