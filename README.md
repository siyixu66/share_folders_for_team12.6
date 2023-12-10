# share_folders_for_team12.6
天体测量及其应用小组作业
selcon.f里的236行换成JPLEPH的绝对路径      NAMFIL='/Users/jcyuan/Desktop/lunar/JPLEPH'
然后编译一下就可生成可执行文件了

test.f90为主程序，可以自行设置开始年份和结束年份
selcon.f里是调用星历及常量的子函数
compute.f是sofa程序库里关于时间转换的函数
