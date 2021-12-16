import os
import time

OPTIONS = ['-lm', '-std=c99', '-Wall', '-O2']
COMPLIE_LIST = [
    'compare',
    'gen_matrix',
    'mult_serial',
    'print'
]

if __name__ == '__main__':
    if os.path.exists('result.stdout'):
        os.remove('result.stdout')
    os.system('mpicc mult_cannon.c -o mult_cannon ' + ' '.join(OPTIONS))
    for i in COMPLIE_LIST:
        print('compiling ' + i)
        os.system('gcc ' + i + '.c -o ' + i + ' ' + ' '.join(OPTIONS))
    n, m, p = input().split()
    os.system('./gen_matrix ' + n + ' ' + m + ' ' + 'matrix_a.stdin')
    os.system('./gen_matrix ' + m + ' ' + p + ' ' + 'matrix_b.stdin')
    if int(n) * int(m) * int(p) < 2048:
        os.system('./print matrix_a.stdin')
        print()
        os.system('./print matrix_b.stdin')
        print()
    os.system('./mult_serial matrix_a.stdin matrix_b.stdin standard.stdout')
    print()
    if int(n) * int(m) * int(p) < 2048:
        os.system('./print standard.stdout')
        print()
    os.system('sbatch run.sh')
    print()
    retry = 0
    while retry < 20 and not os.path.exists('result.stdout'):
        time.sleep(0.5)
        retry += 1
    if retry == 20:
        print('Failed to read result.stdout')
        exit(0)
    if int(n) * int(m) * int(p) < 2048:
        os.system('./print result.stdout')
        print()
    os.system('./compare standard.stdout result.stdout')
    print()
