#!/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/python-2.7.15-tjspyxf3yh64tn4cg34lvdclxruddxaa/bin/python
#

import sys
import os
import subprocess

def cmdlist(str):
    return list(x.strip().replace("'",'') for x in str.split('\n') if x)
env = dict(os.environ)
env['CC'] = '/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-4.8.5/gcc-8.2.0-sxbf4jq6ghmoybsjlpqz2dm2qbbxzfyn/bin/gcc'
env['CMAKE_PREFIX_PATH'] = ":".join(cmdlist("""
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libiconv-1.15-q7w32zbfjms2zfzcllbtovtga4r46u6z
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libxslt-1.1.32-pmogeq5jollvtni2yxzevvhduh4hw3en
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/cmake-3.13.2-lgifhqmxfoaukg3y6ovlxkrfehju6yqy
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openmpi-3.1.3-atscp4jz5a7tt2lig7m2s6tdeb5pgx3t
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/sleef-3.3.1-3swatsby4dy4whgz6ikqykjs7y5qd7ve
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/mpc-1.1.0-k3sl425blztggxpvmjl3ktspwb4vxesr
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/vc-1.3.0-t3mik7lxm3v6tvaoow3tvi7nijcpg4d7
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/sqlite-3.26.0-s5v72k5h3wkzyyf6lodngc367tikautf
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/arblib-2.15.1-vvykvnf37543azw2254ccx66afxda5oe
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openssl-1.1.1-amxovq22o3ukqljhepctty3yvg4xwlzd
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gsl-2.5-zi2ov24dfs5byq5t66nwyxxsw4oiz5me
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/hwloc-1.11.11-6brxyfqbgwgcvsbgbrtfcn4ii2plvf5i
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/bzip2-1.0.6-7djguhm2gowcvvtpweipynfjcj2szvxw
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/catch2-2.5.0-x3alhatojqazushcgx36t53oz7fzmxqv
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/trompeloeil-32-ipzka3enmvb3p4p4mqu5u3xzu4nid775
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/boost-1.69.0-wh2dwx4m5bo6blvx42xjpnszuejnd23i
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/numactl-2.0.12-355ef36k3m3xnqdfmizkys6qa7j6oiss
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/python-2.7.15-tjspyxf3yh64tn4cg34lvdclxruddxaa
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gcc-8.2.0-ofh3irxuvx5t3qq3bhrnrz3mj77ctcss
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/ccache-3.6-f3cx6i5u44bdqjmgvngtuwqbrbsct2p2
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/googletest-1.8.1-tjy6y5jy634pnwahybk4p3dzftcgxu65
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/intel-mkl-2019.1.144-hhawmrvixxoc3gtyathay6j4alvfj6fi
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/doxygen-1.8.15-7qtakpp4w76oats34woylifl6obl66qj
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/zlib-1.2.11-6wdekzyadff6repdb5bxinnzizvnn35t
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/flint-2.5.2-ndhifg24qjt5crxu2nik4mifcsgu4wea
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/xz-5.2.4-fvj6tohizvl45nmccwhd36ppjrtsxrte
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/ncurses-6.1-i655y27j765f243z6i562sryflgenk33
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/intel-tbb-2019.2-oxxqy5lsx6lnx3u3kbcrcip7qz3ymhjn
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gperf-3.0.4-rggielmpgsffz6zp2b7zjdoe7v4fzc2v
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/mpfr-3.1.6-gx2ojhraeqmp5f5xtkqjtp3opahtrftc
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gmp-6.1.2-yhxishziedyfvcfjx7eg6piffz6goeaf
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libgpg-error-1.27-a6z3cplgc7qiar2xjabp2lisq6yemygd
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/hdf5-1.10.4-3quqb62hgt2tgt7xf4wfaaeflc2cjlen
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/erfa-1.4.0-3gm4zhnjhaau7opnv2j5wgxrcraptapw
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/readline-7.0-wcnxwmbecktrveueg3pht7w7vr4c3orn
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libgcrypt-1.8.1-a2nc4vjtt4r5mhlrualfcxbyb2efzgab
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/isl-0.18-idt2taqmn5jr5juzxute6u3qtoxghdho
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/googlebenchmark-1.4.1-2qc3uldtajyhalr7c7sh72zpq3fqm7tr
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libxml2-2.9.8-z2gjacsstfth7p2zx2t2jzdeyk2vmqbo
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libpciaccess-0.13.5-svkyphbquafpjfs4tfsgr3wlmgfknhtu
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gdbm-1.18.1-2ggaajb6mzsqorehy3fzrglcm7t35f3u
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libszip-2.1.1-kg225bg7y5k2y4ouosikcpn7qy5wyzxf
"""))
env['CXX'] = '/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-4.8.5/gcc-8.2.0-sxbf4jq6ghmoybsjlpqz2dm2qbbxzfyn/bin/g++'
env['FC'] = '/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-4.8.5/gcc-8.2.0-sxbf4jq6ghmoybsjlpqz2dm2qbbxzfyn/bin/gfortran'
env['PATH'] = ":".join(cmdlist("""
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/python-2.7.15-tjspyxf3yh64tn4cg34lvdclxruddxaa/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gsl-2.5-zi2ov24dfs5byq5t66nwyxxsw4oiz5me/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openssl-1.1.1-amxovq22o3ukqljhepctty3yvg4xwlzd/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/cmake-3.13.2-lgifhqmxfoaukg3y6ovlxkrfehju6yqy/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gcc-8.2.0-ofh3irxuvx5t3qq3bhrnrz3mj77ctcss/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/hdf5-1.10.4-3quqb62hgt2tgt7xf4wfaaeflc2cjlen/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/python-2.7.15-tjspyxf3yh64tn4cg34lvdclxruddxaa/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/ccache-3.6-f3cx6i5u44bdqjmgvngtuwqbrbsct2p2/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/graphviz-2.40.1-wujmug3qnpb47hwcdr4x3q3oby4p6kcs/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/doxygen-1.8.15-7qtakpp4w76oats34woylifl6obl66qj/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/intel-mkl-2019.1.144-hhawmrvixxoc3gtyathay6j4alvfj6fi/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gsl-2.5-zi2ov24dfs5byq5t66nwyxxsw4oiz5me/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openssl-1.1.1-amxovq22o3ukqljhepctty3yvg4xwlzd/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/cmake-3.13.2-lgifhqmxfoaukg3y6ovlxkrfehju6yqy/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gcc-8.2.0-ofh3irxuvx5t3qq3bhrnrz3mj77ctcss/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/hdf5-1.10.4-3quqb62hgt2tgt7xf4wfaaeflc2cjlen/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/python-2.7.15-tjspyxf3yh64tn4cg34lvdclxruddxaa/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/ccache-3.6-f3cx6i5u44bdqjmgvngtuwqbrbsct2p2/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/doxygen-1.8.15-7qtakpp4w76oats34woylifl6obl66qj/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/intel-mkl-2019.1.144-hhawmrvixxoc3gtyathay6j4alvfj6fi/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/pricer-development-4s55oak2p3cvmdudcutsq36ag4s4dx3m/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/intel-mkl-2019.1.144-hhawmrvixxoc3gtyathay6j4alvfj6fi/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/hdf5-1.10.4-3quqb62hgt2tgt7xf4wfaaeflc2cjlen/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gsl-2.5-zi2ov24dfs5byq5t66nwyxxsw4oiz5me/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gcc-8.2.0-ofh3irxuvx5t3qq3bhrnrz3mj77ctcss/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/doxygen-1.8.15-7qtakpp4w76oats34woylifl6obl66qj/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/graphviz-2.40.1-wujmug3qnpb47hwcdr4x3q3oby4p6kcs/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/cmake-3.13.2-lgifhqmxfoaukg3y6ovlxkrfehju6yqy/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/ccache-3.6-f3cx6i5u44bdqjmgvngtuwqbrbsct2p2/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libxslt-1.1.32-pmogeq5jollvtni2yxzevvhduh4hw3en/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libgcrypt-1.8.1-a2nc4vjtt4r5mhlrualfcxbyb2efzgab/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libgpg-error-1.27-a6z3cplgc7qiar2xjabp2lisq6yemygd/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gperf-3.0.4-rggielmpgsffz6zp2b7zjdoe7v4fzc2v/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/python-2.7.15-tjspyxf3yh64tn4cg34lvdclxruddxaa/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/sqlite-3.26.0-s5v72k5h3wkzyyf6lodngc367tikautf/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openssl-1.1.1-amxovq22o3ukqljhepctty3yvg4xwlzd/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gdbm-1.18.1-2ggaajb6mzsqorehy3fzrglcm7t35f3u/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/readline-7.0-wcnxwmbecktrveueg3pht7w7vr4c3orn/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/ncurses-6.1-i655y27j765f243z6i562sryflgenk33/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openmpi-3.1.3-atscp4jz5a7tt2lig7m2s6tdeb5pgx3t/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/hwloc-1.11.11-6brxyfqbgwgcvsbgbrtfcn4ii2plvf5i/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/numactl-2.0.12-355ef36k3m3xnqdfmizkys6qa7j6oiss/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libxml2-2.9.8-z2gjacsstfth7p2zx2t2jzdeyk2vmqbo/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/xz-5.2.4-fvj6tohizvl45nmccwhd36ppjrtsxrte/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libiconv-1.15-q7w32zbfjms2zfzcllbtovtga4r46u6z/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/bzip2-1.0.6-7djguhm2gowcvvtpweipynfjcj2szvxw/bin
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-4.8.5/environment-modules-3.2.10-k5d3cqqywkaxlplgqilode6hcrwmzjix/Modules/bin
    /home/VortexUser/spack/bin
    /home/VortexUser/.opam/system/bin
    /home/VortexUser/anaconda3/bin
    /usr/local/bin
    /usr/local/sbin
    /usr/bin
    /usr/sbin
    /bin
    /sbin
    /home/VortexUser/.local/bin
    /home/VortexUser/bin
    /home/VortexUser/ocamlbrew/ocaml-4.05.0/bin
"""))
env['SPACK_TRANSITIVE_INCLUDE_PATH'] = ";".join(cmdlist("""
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libiconv-1.15-q7w32zbfjms2zfzcllbtovtga4r46u6z/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libxslt-1.1.32-pmogeq5jollvtni2yxzevvhduh4hw3en/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/cmake-3.13.2-lgifhqmxfoaukg3y6ovlxkrfehju6yqy/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openmpi-3.1.3-atscp4jz5a7tt2lig7m2s6tdeb5pgx3t/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/sleef-3.3.1-3swatsby4dy4whgz6ikqykjs7y5qd7ve/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/mpc-1.1.0-k3sl425blztggxpvmjl3ktspwb4vxesr/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/vc-1.3.0-t3mik7lxm3v6tvaoow3tvi7nijcpg4d7/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/sqlite-3.26.0-s5v72k5h3wkzyyf6lodngc367tikautf/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/arblib-2.15.1-vvykvnf37543azw2254ccx66afxda5oe/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openssl-1.1.1-amxovq22o3ukqljhepctty3yvg4xwlzd/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gsl-2.5-zi2ov24dfs5byq5t66nwyxxsw4oiz5me/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/hwloc-1.11.11-6brxyfqbgwgcvsbgbrtfcn4ii2plvf5i/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/bzip2-1.0.6-7djguhm2gowcvvtpweipynfjcj2szvxw/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/catch2-2.5.0-x3alhatojqazushcgx36t53oz7fzmxqv/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/trompeloeil-32-ipzka3enmvb3p4p4mqu5u3xzu4nid775/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/boost-1.69.0-wh2dwx4m5bo6blvx42xjpnszuejnd23i/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/numactl-2.0.12-355ef36k3m3xnqdfmizkys6qa7j6oiss/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/python-2.7.15-tjspyxf3yh64tn4cg34lvdclxruddxaa/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gcc-8.2.0-ofh3irxuvx5t3qq3bhrnrz3mj77ctcss/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/ccache-3.6-f3cx6i5u44bdqjmgvngtuwqbrbsct2p2/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/googletest-1.8.1-tjy6y5jy634pnwahybk4p3dzftcgxu65/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/intel-mkl-2019.1.144-hhawmrvixxoc3gtyathay6j4alvfj6fi/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/doxygen-1.8.15-7qtakpp4w76oats34woylifl6obl66qj/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/zlib-1.2.11-6wdekzyadff6repdb5bxinnzizvnn35t/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/flint-2.5.2-ndhifg24qjt5crxu2nik4mifcsgu4wea/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/xz-5.2.4-fvj6tohizvl45nmccwhd36ppjrtsxrte/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/ncurses-6.1-i655y27j765f243z6i562sryflgenk33/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/intel-tbb-2019.2-oxxqy5lsx6lnx3u3kbcrcip7qz3ymhjn/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gperf-3.0.4-rggielmpgsffz6zp2b7zjdoe7v4fzc2v/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/mpfr-3.1.6-gx2ojhraeqmp5f5xtkqjtp3opahtrftc/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gmp-6.1.2-yhxishziedyfvcfjx7eg6piffz6goeaf/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libgpg-error-1.27-a6z3cplgc7qiar2xjabp2lisq6yemygd/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/hdf5-1.10.4-3quqb62hgt2tgt7xf4wfaaeflc2cjlen/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/erfa-1.4.0-3gm4zhnjhaau7opnv2j5wgxrcraptapw/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/readline-7.0-wcnxwmbecktrveueg3pht7w7vr4c3orn/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libgcrypt-1.8.1-a2nc4vjtt4r5mhlrualfcxbyb2efzgab/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/isl-0.18-idt2taqmn5jr5juzxute6u3qtoxghdho/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/googlebenchmark-1.4.1-2qc3uldtajyhalr7c7sh72zpq3fqm7tr/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libxml2-2.9.8-z2gjacsstfth7p2zx2t2jzdeyk2vmqbo/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libpciaccess-0.13.5-svkyphbquafpjfs4tfsgr3wlmgfknhtu/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gdbm-1.18.1-2ggaajb6mzsqorehy3fzrglcm7t35f3u/include
    /home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libszip-2.1.1-kg225bg7y5k2y4ouosikcpn7qy5wyzxf/include
"""))

cmd = cmdlist("""
/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/cmake-3.13.2-lgifhqmxfoaukg3y6ovlxkrfehju6yqy/bin/cmake
    -G
    Unix Makefiles
    -DCMAKE_INSTALL_PREFIX:PATH=/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/pricer-development-4s55oak2p3cvmdudcutsq36ag4s4dx3m
    -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=FALSE
    -DCMAKE_INSTALL_RPATH:STRING=/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/pricer-development-4s55oak2p3cvmdudcutsq36ag4s4dx3m/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/pricer-development-4s55oak2p3cvmdudcutsq36ag4s4dx3m/lib64;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/arblib-2.15.1-vvykvnf37543azw2254ccx66afxda5oe/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/flint-2.5.2-ndhifg24qjt5crxu2nik4mifcsgu4wea/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gmp-6.1.2-yhxishziedyfvcfjx7eg6piffz6goeaf/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/mpfr-3.1.6-gx2ojhraeqmp5f5xtkqjtp3opahtrftc/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/boost-1.69.0-wh2dwx4m5bo6blvx42xjpnszuejnd23i/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/bzip2-1.0.6-7djguhm2gowcvvtpweipynfjcj2szvxw/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openmpi-3.1.3-atscp4jz5a7tt2lig7m2s6tdeb5pgx3t/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/hwloc-1.11.11-6brxyfqbgwgcvsbgbrtfcn4ii2plvf5i/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libpciaccess-0.13.5-svkyphbquafpjfs4tfsgr3wlmgfknhtu/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libxml2-2.9.8-z2gjacsstfth7p2zx2t2jzdeyk2vmqbo/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libiconv-1.15-q7w32zbfjms2zfzcllbtovtga4r46u6z/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/xz-5.2.4-fvj6tohizvl45nmccwhd36ppjrtsxrte/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/zlib-1.2.11-6wdekzyadff6repdb5bxinnzizvnn35t/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/numactl-2.0.12-355ef36k3m3xnqdfmizkys6qa7j6oiss/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/python-2.7.15-tjspyxf3yh64tn4cg34lvdclxruddxaa/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gdbm-1.18.1-2ggaajb6mzsqorehy3fzrglcm7t35f3u/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/readline-7.0-wcnxwmbecktrveueg3pht7w7vr4c3orn/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/ncurses-6.1-i655y27j765f243z6i562sryflgenk33/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openssl-1.1.1-amxovq22o3ukqljhepctty3yvg4xwlzd/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/sqlite-3.26.0-s5v72k5h3wkzyyf6lodngc367tikautf/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libxslt-1.1.32-pmogeq5jollvtni2yxzevvhduh4hw3en/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libgcrypt-1.8.1-a2nc4vjtt4r5mhlrualfcxbyb2efzgab/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libgpg-error-1.27-a6z3cplgc7qiar2xjabp2lisq6yemygd/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/erfa-1.4.0-3gm4zhnjhaau7opnv2j5wgxrcraptapw/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gcc-8.2.0-ofh3irxuvx5t3qq3bhrnrz3mj77ctcss/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/isl-0.18-idt2taqmn5jr5juzxute6u3qtoxghdho/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/mpc-1.1.0-k3sl425blztggxpvmjl3ktspwb4vxesr/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/googlebenchmark-1.4.1-2qc3uldtajyhalr7c7sh72zpq3fqm7tr/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gsl-2.5-zi2ov24dfs5byq5t66nwyxxsw4oiz5me/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/hdf5-1.10.4-3quqb62hgt2tgt7xf4wfaaeflc2cjlen/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/libszip-2.1.1-kg225bg7y5k2y4ouosikcpn7qy5wyzxf/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/intel-mkl-2019.1.144-hhawmrvixxoc3gtyathay6j4alvfj6fi/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/intel-tbb-2019.2-oxxqy5lsx6lnx3u3kbcrcip7qz3ymhjn/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/sleef-3.3.1-3swatsby4dy4whgz6ikqykjs7y5qd7ve/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/trompeloeil-32-ipzka3enmvb3p4p4mqu5u3xzu4nid775/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/vc-1.3.0-t3mik7lxm3v6tvaoow3tvi7nijcpg4d7/lib;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/catch2-2.5.0-x3alhatojqazushcgx36t53oz7fzmxqv/lib64;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gcc-8.2.0-ofh3irxuvx5t3qq3bhrnrz3mj77ctcss/lib64;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/googletest-1.8.1-tjy6y5jy634pnwahybk4p3dzftcgxu65/lib64
    -DCMAKE_PREFIX_PATH:STRING=/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gcc-8.2.0-ofh3irxuvx5t3qq3bhrnrz3mj77ctcss;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/doxygen-1.8.15-7qtakpp4w76oats34woylifl6obl66qj;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/cmake-3.13.2-lgifhqmxfoaukg3y6ovlxkrfehju6yqy;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/googlebenchmark-1.4.1-2qc3uldtajyhalr7c7sh72zpq3fqm7tr;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/hdf5-1.10.4-3quqb62hgt2tgt7xf4wfaaeflc2cjlen;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/python-2.7.15-tjspyxf3yh64tn4cg34lvdclxruddxaa;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/intel-tbb-2019.2-oxxqy5lsx6lnx3u3kbcrcip7qz3ymhjn;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/catch2-2.5.0-x3alhatojqazushcgx36t53oz7fzmxqv;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/sleef-3.3.1-3swatsby4dy4whgz6ikqykjs7y5qd7ve;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/openssl-1.1.1-amxovq22o3ukqljhepctty3yvg4xwlzd;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/trompeloeil-32-ipzka3enmvb3p4p4mqu5u3xzu4nid775;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/vc-1.3.0-t3mik7lxm3v6tvaoow3tvi7nijcpg4d7;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/ccache-3.6-f3cx6i5u44bdqjmgvngtuwqbrbsct2p2;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/erfa-1.4.0-3gm4zhnjhaau7opnv2j5wgxrcraptapw;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/arblib-2.15.1-vvykvnf37543azw2254ccx66afxda5oe;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/intel-mkl-2019.1.144-hhawmrvixxoc3gtyathay6j4alvfj6fi;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gmp-6.1.2-yhxishziedyfvcfjx7eg6piffz6goeaf;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/boost-1.69.0-wh2dwx4m5bo6blvx42xjpnszuejnd23i;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/gsl-2.5-zi2ov24dfs5byq5t66nwyxxsw4oiz5me;/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/mpfr-3.1.6-gx2ojhraeqmp5f5xtkqjtp3opahtrftc
    BOOST_DIR=/home/VortexUser/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/boost-1.69.0-wh2dwx4m5bo6blvx42xjpnszuejnd23i
""") + sys.argv[1:]

proc = subprocess.Popen(cmd, env=env)
proc.wait()
