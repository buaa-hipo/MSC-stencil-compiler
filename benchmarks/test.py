
from os import system as cmd
import os
import time

#stencils_2d = ["2d9pt_star", "2d9pt_box", "2d121pt_box", "2d169pt_box"]
stencils_2d = []
stencils_3d = ["3d7pt_star", "3d13pt_star", "3d25pt_star", "3d31pt_star"]
#stencils_3d = ["3d25pt_star", "3d31pt_star"]
stencils_3d = ["3d13pt_star"]
stencils = stencils_3d + stencils_2d

tile_size_dict = {"2d9pt_star":[32,64], "2d9pt_box":[32,64], 
                  "2d121pt_box":[16,32], "2d169pt_box":[16,32],
                  "3d7pt_star":[2,8,64], "3d13pt_star":[2,8,64], 
                  "3d25pt_star":[2,4,32], "3d31pt_star":[2,4,32] }

def get_file(substr, lists):
    for name in lists:
        if substr in name:
            return name
    return "not found requested file"

def check_file(substr, lists):
    for name in lists:
        if substr in name:
            return True
    return False

def get_mul(lists):
    res=1
    for ele in lists:
        res=res*ele
    return res

def test(bench_kind, target, num_threads, cases, config, path_out, num_mpi, out_flag=""):
    cmd("make clear")
    #config = bench_kind+" "+target+" "+str(num_threads)+" "+config
    cmd("./benchmark "+config)
    files = os.listdir("./")
    compile_file = get_file("compile", files)
    #
    run_file =""
    if(target=="feiteng"):
        run_file = get_file("batch", files)
    elif(target=="sunway"):
        run_file = get_file("run", files)
    elif(target=="x86"):
        run_file = get_file("run", files)
    cmd("sh "+compile_file+ " && "+"sh "+ run_file)
    #
    part_out_file=""
    if(target=="feiteng"):
        part_out_file = "slurm"
    elif(target=="sunway"):
        part_out_file = ".out"
    elif(target=="x86"):
        part_out_file = ".x86.out"
    #
    while True:
        files = os.listdir("./")
        if(check_file(part_out_file,files)) :
            print("find "+part_out_file)
            #time.sleep(60)
            files = os.listdir("./")
            out_file = get_file(part_out_file, files)
            f = open(out_file, 'r')
            content = f.read()
            #print(content)
            if "mpi" in content:
                print("find mpi")
                break
            else:
                print("no find mpi")
                time.sleep(1)
                f.close()
                continue
            #
        print("not find "+part_out_file)
        time.sleep(1)
    #
    files = os.listdir("./")
    out_file = get_file(part_out_file, files)
    print(path_out)
    if not os.path.exists(path_out):
        print(path_out+" not exist!")
        os.makedirs(path_out)
    print("mv "+out_file+" "+path_out+out_file+"."+str(num_mpi)+out_flag)
    cmd("mv "+out_file+" "+path_out+out_file+"."+str(num_mpi)+out_flag)


def test_scalability(bench_kind, flag_scalability, target, num_threads, shape_grid_3d, shape_sub_grid_3d, tile_size_xyz, shape_grid_2d, shape_sub_grid_2d, tile_size_xy, out_flag=""):
    #
    for cases in stencils:
        print("cases: "+str(cases))
        path_out = "./out/"+target+"/"+flag_scalability+"/"+str(cases)+"/"
        print(path_out)
        if(cases in stencils_3d):
            for i in range(len(shape_grid_3d)):
                config = str(cases)+" "
                config += target+" "
                config += str(num_threads[i][0])+" "
                config += str(shape_grid_3d[i][0])+" "
                config += str(shape_grid_3d[i][1])+" "
                config += str(shape_grid_3d[i][2])+" "
                config += str(shape_sub_grid_3d[i][0])+" "
                config += str(shape_sub_grid_3d[i][1])+" "
                config += str(shape_sub_grid_3d[i][2])+" "
                if not (target == "sunway"):
                    config += str(tile_size_xyz[i][0])+" "
                    config += str(tile_size_xyz[i][1])+" "
                    config += str(tile_size_xyz[i][2])+" "
                else:
                    config += str(tile_size_dict[cases][0])+" "
                    config += str(tile_size_dict[cases][1])+" "
                    config += str(tile_size_dict[cases][2])+" "
                config += str(1)+" " # use_schedule
                config += str(0)+" " # is_assigned
                #
                print("config: "+str(config))
                #continue
                test(bench_kind, target, num_threads[i][0], cases, config, path_out, get_mul(shape_grid_3d[i]), out_flag)
        if(cases in stencils_2d):
            for i in range(len(shape_grid_2d)):
                config = str(cases)+" "
                config += target+" "
                config += str(num_threads[i][0])+" "
                config += str(shape_grid_2d[i][0])+" "
                config += str(shape_grid_2d[i][1])+" "
                config += str(shape_sub_grid_2d[i][0])+" "
                config += str(shape_sub_grid_2d[i][1])+" "
                if not (target == "sunway"):
                    config += str(tile_size_xy[i][0])+" "
                    config += str(tile_size_xy[i][1])+" "
                else:
                    config += str(tile_size_dict[cases][0])+" "
                    config += str(tile_size_dict[cases][1])+" "
                config += str(1)+" " # use_schedule
                config += str(0)+" " # is_assigned
                #
                #config = flag_scalability + " " + config
                print("config: "+str(config))
                #continue
                test(bench_kind, target, num_threads[i][0], cases, config, path_out, get_mul(shape_grid_2d[i]), out_flag)

def test_strong_scalability(target):
    if (target == "feiteng"):
        shape_grid_3d=[[4,4,2], [4,4,4], [4,8,4], [8,8,4]]
        shape_sub_grid_3d = [[256,256,256], [256,256,128], [256,128,128], [128,128,128]]
        tile_size_xyz = [[2,16,128], [2,16,128], [2,16,128], [2,16,128]]
        shape_grid_2d=[[8,4], [8,8], [16,8], [16,16]]
        shape_sub_grid_2d = [[4096,4096], [4096,2048], [2048,2048], [2048,1024]]
        tile_size_xy = [[4,1024], [4,1024], [4,1024], [4,1024]]
        num_threads = [[32], [32], [32], [32]]
        #
        test_scalability("scalability", "strong_scalability", target, num_threads, shape_grid_3d, shape_sub_grid_3d,
                tile_size_xyz, shape_grid_2d, shape_sub_grid_2d, tile_size_xy)
    elif (target == "sunway"):
        shape_grid_3d=[[8,4,4], [8,8,4], [8,8,8], [16,8,8]]
        shape_sub_grid_3d = [[256,512,512], [256,256,512], [256,256,256], [128,256,256]]
        tile_size_xyz = None
        shape_grid_2d=[[8,16], [16,16], [16,32], [32,32]]
        shape_sub_grid_2d = [[8192,4096], [4096,4096], [4096,2048], [2048,2048]]
        tile_size_xy = None
        num_threads = [[64], [64], [64], [64]]
        #
        test_scalability("scalability", "strong_scalability", target, num_threads, shape_grid_3d, shape_sub_grid_3d,\
                tile_size_xyz, shape_grid_2d, shape_sub_grid_2d, tile_size_xy)
    else:
        print("error ! not correct target!")
        exit()

def test_weak_scalability(target):
    if (target == "feiteng"):
        shape_grid_3d=[[4,4,2], [4,4,4], [4,8,4], [8,8,4]]
        shape_sub_grid_3d = [[256,256,256], [256,256,256], [256,256,256], [256,256,256]]
        tile_size_xyz = [[2,8,256], [2,8,256], [2,8,256], [2,8,256]]
        shape_grid_2d=[[8,4], [8,8], [16,8], [16,16]]
        shape_sub_grid_2d = [[4096,4096], [4096,4096], [4096,4096], [4096,4096]]
        tile_size_xy = [[2,2048], [2,2048], [2,2048], [2,2048]]
        num_threads = [[32], [32], [32], [32]]
        #
        test_scalability("scalability", "weak_scalability", target, num_threads, shape_grid_3d, shape_sub_grid_3d,\
                tile_size_xyz, shape_grid_2d, shape_sub_grid_2d, tile_size_xy)
    elif (target == "sunway"):
        shape_grid_3d=[[8,4,4], [8,8,4], [8,8,8], [16,8,8]]
        shape_sub_grid_3d = [[256,256,256] for i in range(4)]
        tile_size_xyz = None
        shape_grid_2d=[[16,8], [16,16], [32,16], [32,32]]
        shape_sub_grid_2d = [[4096,4096] for i in range(4)]
        tile_size_xy = None
        num_threads = [[64], [64], [64], [64]]
        #
        test_scalability("scalability", "weak_scalability", target, num_threads, shape_grid_3d, shape_sub_grid_3d,\
                tile_size_xyz, shape_grid_2d, shape_sub_grid_2d, tile_size_xy)
    else:
        print("error ! not correct target!")
        exit()

def test_x86_vs_physis(target):
    if (target == "x86"):
        shape_grid_3d=[[2,2,7], [1,2,7], [1,1,7]]
        shape_sub_grid_3d = [[256,256,256], [512,256,256], [512,512,256]]
        tile_size_xyz = [[2,8,256], [2,8,256], [2,8,256]]
        #
        shape_grid_2d=[[4,7], [2,7], [1,7]]
        shape_sub_grid_2d = [[4096,4096], [8192,4096], [16384,4096]]
        tile_size_xy = [[2,2048], [2,2048], [2,2048]]
        #
        num_threads = [[1], [2], [4]]
        test_scalability("scalability", "dsl_vs_physis", target, num_threads, shape_grid_3d, shape_sub_grid_3d,\
                tile_size_xyz, shape_grid_2d, shape_sub_grid_2d, tile_size_xy)
    else:
        print("error ! not correct target in test x86 vs physis!")
        exit()

def bench(bench_kind, flag_scalability, target, num_threads, shape_grid_3d, shape_sub_grid_3d, tile_size_xyz, shape_grid_2d, shape_sub_grid_2d, tile_size_xy, use_schedule, use_tile, use_reorder, use_threads, is_assigned, out_flag):
    for cases in stencils:
        print("cases: "+str(cases))
        path_out = "./out/"+target+"/"+flag_scalability+"/"+str(cases)+"/"
        print(path_out)
        if(cases in stencils_3d):
            for i in range(len(shape_grid_3d)):
                config = str(cases)+" "
                config += target+" "
                config += str(num_threads[i][0])+" "
                config += str(shape_grid_3d[i][0])+" "
                config += str(shape_grid_3d[i][1])+" "
                config += str(shape_grid_3d[i][2])+" "
                config += str(shape_sub_grid_3d[i][0])+" "
                config += str(shape_sub_grid_3d[i][1])+" "
                config += str(shape_sub_grid_3d[i][2])+" "
                if not (target == "sunway"):
                    config += str(tile_size_xyz[i][0])+" "
                    config += str(tile_size_xyz[i][1])+" "
                    config += str(tile_size_xyz[i][2])+" "
                else:
                    config += str(tile_size_dict[cases][0])+" "
                    config += str(tile_size_dict[cases][1])+" "
                    config += str(tile_size_dict[cases][2])+" "
                config += str(use_schedule)+" " # use_schedule
                config += str(is_assigned)+" " # is_assigned
                #
                print("config: "+str(config))
                #continue
                #test(bench_kind, target, num_threads[i][0], cases, config, path_out, get_mul(shape_grid_3d[i]), out_flag)
        if(cases in stencils_2d):
            for i in range(len(shape_grid_2d)):
                config = str(cases)+" "
                config += target+" "
                config += str(num_threads[i][0])+" "
                config += str(shape_grid_2d[i][0])+" "
                config += str(shape_grid_2d[i][1])+" "
                config += str(shape_sub_grid_2d[i][0])+" "
                config += str(shape_sub_grid_2d[i][1])+" "
                if not (target == "sunway"):
                    config += str(tile_size_xy[i][0])+" "
                    config += str(tile_size_xy[i][1])+" "
                else:
                    config += str(tile_size_dict[cases][0])+" "
                    config += str(tile_size_dict[cases][1])+" "
                config += str(use_schedule)+" " # use_schedule
                config += str(is_assigned)+" " # is_assigned
                #
                print("config: "+str(config))
                #continue
                #test(bench_kind, target, num_threads[i][0], cases, config, path_out, get_mul(shape_grid_2d[i]), out_flag)

def test_x86_patus(target):
    if (target == "x86"):
        shape_grid_3d=[[1,1,1]]
        shape_sub_grid_3d = [[256,256,256]]
        tile_size_xyz = [[2,8,256]]
        shape_grid_2d=[[1,1]]
        shape_sub_grid_2d = [[4096,4096]]
        tile_size_xy = [[2,2048]]
        num_threads = [[28]]
        #
        use_schedule=1 
        use_tile=1
        use_reorder=1 
        use_threads=0 
        is_assigned=0
        #
        assert(not (use_threads == 0) )
        bench("bench_schedule", "dsl_vs_patus", target, num_threads, shape_grid_3d, shape_sub_grid_3d,\
            tile_size_xyz, shape_grid_2d, shape_sub_grid_2d, tile_size_xy,\
            use_schedule, use_tile, use_reorder, use_threads, is_assigned, ".schedule_tile_reorder"+"."+"f64")


def bench_sunway(bench_kind, flag_scalability, target, num_threads, shape_grid_3d, shape_sub_grid_3d, tile_size_xyz, shape_grid_2d, shape_sub_grid_2d, tile_size_xy, use_schedule, is_assigned, out_flag):
    for cases in stencils:
        print("cases: "+str(cases))
        path_out = "./out/"+target+"/"+flag_scalability+"/"+str(cases)+"/"
        print(path_out)
        if(cases in stencils_3d):
            for i in range(len(shape_grid_3d)):
                config = str(cases)+" "
                config += target+" "
                config += str(num_threads[i][0])+" "
                config += str(shape_grid_3d[i][0])+" "
                config += str(shape_grid_3d[i][1])+" "
                config += str(shape_grid_3d[i][2])+" "
                config += str(shape_sub_grid_3d[i][0])+" "
                config += str(shape_sub_grid_3d[i][1])+" "
                config += str(shape_sub_grid_3d[i][2])+" "
                config += str(tile_size_dict[cases][0])+" "
                config += str(tile_size_dict[cases][1])+" "
                config += str(tile_size_dict[cases][2])+" "
                config += str(use_schedule)+" " # use_schedule
                config += str(is_assigned)+" " # is_assigned

                print("config: "+str(config))
                #continue
                test(bench_kind, target, num_threads[i][0], cases, config, path_out, get_mul(shape_grid_3d[i]), out_flag)
        if(cases in stencils_2d):
            for i in range(len(shape_grid_2d)):
                config = str(cases)+" "
                config += target+" "
                config += str(num_threads[i][0])+" "
                config += str(shape_grid_2d[i][0])+" "
                config += str(shape_grid_2d[i][1])+" "
                config += str(shape_sub_grid_2d[i][0])+" "
                config += str(shape_sub_grid_2d[i][1])+" "
                config += str(tile_size_dict[cases][0])+" "
                config += str(tile_size_dict[cases][1])+" "
                config += str(use_schedule)+" " # use_schedule
                config += str(is_assigned)+" " # is_assigned

                print("config: "+str(config))
                #continue
                test(bench_kind, target, num_threads[i][0], cases, config, path_out, get_mul(shape_grid_2d[i]), out_flag)

def test_sunway(target, bench_flag="bench_sunway"):
    shape_grid_3d=[[1,1,1]]
    shape_sub_grid_3d = [[256,256,256]]
    tile_size_xyz = None
    shape_grid_2d=[[1,1]]
    shape_sub_grid_2d = [[4096,4096]]
    tile_size_xy = None
    num_threads = [[64]]
    # 
    if (target == "sunway" and bench_flag == "bench_sunway"):
        #
        use_schedule=1
        is_assigned=0
        #
        bench_sunway("bench_sunway", "dsl_vs_openacc", target, num_threads, shape_grid_3d, shape_sub_grid_3d,\
                tile_size_xyz, shape_grid_2d, shape_sub_grid_2d, tile_size_xy, use_schedule, is_assigned, ".not_assigned")


def main():
    #cmd("g++ benchmark.cc -std=c++11 -o benchmark")
    target="sunway"
    #test_sunway(target)
    test_weak_scalability(target)
    test_strong_scalability(target)
    target="x86"
    #test_x86_patus(target)
    #test_x86_vs_physis(target)
    target="feiteng"
    #

if __name__ == "__main__":
    main()

'''
##
shape_grid_3d32mpi  = [4,4,2]
shape_grid_3d64mpi  = [4,4,4]
shape_grid_3d128mpi = [4,8,4]
shape_grid_3d256mpi = [8,8,4]
shape_grid_3d=[[4,4,2], [4,4,4], [4,8,4], [8,8,4]]
##
shape_sub_grid_3d32mpi  = [256,256,256]
shape_sub_grid_3d64mpi  = [256,256,128]
shape_sub_grid_3d128mpi = [256,128,128]
shape_sub_grid_3d256mpi = [128,128,128]
shape_sub_grid_3d = [[256,256,256], [256,256,128], [256,128,128], [128,128,128]]
##
tile_size_xyz_3d32mpi  = [2,16,128]
tile_size_xyz_3d64mpi  = [2,16,128]
tile_size_xyz_3d128mpi = [2,16,128]
tile_size_xyz_3d256mpi = [2,16,128]
tile_size_xyz = [[2,16,128], [2,16,128], [2,16,128], [2,16,128]]
##
shape_grid_2d=[[8,4], [8,8], [16,8], [16,16]]
##
shape_sub_grid_2d = [[4096,4096], [4096,2048], [2048,2048], [2048,1024]]
##
tile_size_xy = [[4,1024], [4,1024], [4,1024], [4,1024]]
'''

'''
def test_x86_vs_physis(target):
    if (target == "x86"):
        shape_grid_3d=[[2,2,7], [1,2,7], [1,1,7], [2,2,1]]
        shape_sub_grid_3d = [[256,256,256], [512,256,256], [512,512,256], [256,256,1792]]
        tile_size_xyz = [[2,8,256], [2,8,256], [2,8,256], [2,8,256]]
        #
        shape_grid_2d=[[4,7], [2,7], [1,7], [4,1]]
        shape_sub_grid_2d = [[4096,4096], [8192,4096], [16384,4096], [4096,28672]]
        tile_size_xy = [[2,2048], [2,2048], [2,2048], [2,2048]]
        #
        num_threads = [[1], [2], [4], [7]]
        test_scalability("dsl_vs_physis", target, num_threads, shape_grid_3d, shape_sub_grid_3d,\
                tile_size_xyz, shape_grid_2d, shape_sub_grid_2d, tile_size_xy)
    else:
        print("error ! not correct target in test x86 vs physis!")
        exit()
'''


