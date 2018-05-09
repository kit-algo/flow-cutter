#!/usr/bin/env python2
import sys, os, re, shutil, subprocess

################################################################################
## Simple Configuration
##

# Use a color terminal?
use_colors = True

# What compiler to use?
compiler = "g++"
#compiler = "clang++"

# Use C++11 ?
use_cpp11 = True

# Be verbose ?
use_verbose = ("--verbose" in sys.argv)

# Clean up after being finished?
use_clean = ("--clean" in sys.argv)

# Ignore warning
ignore_warnings = ("--ignore-warnings" in sys.argv)

# Get rid of GPL dependencies
no_gpl = ("--no-gpl" in sys.argv)

# Clean up after being finished?
show_header_scanning = ("--show-header-scanning" in sys.argv)

# Commands to add to each file
# "-D_GLIBCXX_DEBUG"
# "-msse4a", "-DBE_COMPUTE2_COMPATIBLE" "-O3", "-DNDEBUG"
compiler_settings = ["-O3", "-DNDEBUG"]
if ignore_warnings:
	compiler_settings += ["-w"]
else:
	compiler_settings += ["-Wall", "-Wdisabled-optimization"]

if no_gpl:
	compiler_settings += ["-DNO_GPL"]

linker_settings = []

source_extensions = [".cpp", ".cxx"]
header_extensions = [".h", ".hpp", ".hxx"]

################################################################################
## Advanced Configuration
##

if use_cpp11:
	compiler_settings.append("-std=c++0x")
	linker_settings.append("-std=c++0x")

def make_object_file_name(whole_path):
	path, file_name = os.path.split(whole_path)
	base_name, extension = os.path.splitext(file_name)
	return os.path.join(path, "."+base_name+".o")

if use_colors:
	color_command = '\033[92m'
	color_error = '\033[91m'
	color_end = '\033[0m'
else:
	color_command = ''
	color_error = ''
	color_end = ''

################################################################################
## Basic commands and source file information extraction
##

def invisible_cmd(command):
	p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = p.communicate()
	if p.returncode != 0 or stderr != "":
		print(color_error+"Command Failed:" + color_end + "\n" + stderr)
		sys.exit(1)
	return stdout

def cmd(command):
	print(color_command+" ".join(command)+color_end)
	p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = p.communicate()
	if p.returncode != 0 or stderr != "":
		print(color_error+"Command Failed:" + color_end + "\n" + stderr)
		sys.exit(1)
	return stdout

def find_all_includes(source_file):
	c = invisible_cmd
	if show_header_scanning:
		c = cmd
	return list(x for x in re.findall("(\\S+)",
		c([compiler, "-E", "-M"] + compiler_settings + [source_file]))
		if x != "\\" and not os.path.splitext(x)[1] in source_extensions + [".o:"]
	)


def run_compiler(source_file, object_file, flags):
	cmd([compiler, "-c"] + compiler_settings + flags + [source_file, "-o", object_file])

def run_linker(object_files, executable_file, flags):
	cmd([compiler] + object_files + ["-o", executable_file] + linker_settings + flags)

def find_flags(code, flag_name):
	l = (re.split('\s+', x) for x in re.findall(re.escape(flag_name)+'(.*)', code))
	return list(x for y in l for x in y if x!="")

def extract_code_flags(file_name):
	with open(file_name,"r") as f:
		code = f.read()
		return (
			find_flags(code, "compile with"),
			find_flags(code, "compile related with"),
			find_flags(code, "compile all with"),
			find_flags(code, "link with")
		)

def get_last_modified(file_path):
	if os.path.exists(file_path):
		return os.path.getmtime(file_path)
	else:
		return 0

################################################################################
## Find all source files
##
main_file_list = []

def gen_file_lists():
	global source_file_list
	global main_file_list
	global header_file_list

	header_file_list = []
	source_file_list = []

	for root, dirs, files in os.walk("."):
		for file_name in files:
			if file_name[0]!= ".":
				ext = os.path.splitext(file_name)[1]
				if(ext in source_extensions):
					source_file_list.append(os.path.relpath(os.path.join(root, file_name), "."))
				if(ext in header_extensions):
					header_file_list.append(os.path.relpath(os.path.join(root, file_name), "."))


	main_file_list = list(x for x in source_file_list
		if
			re.search(r"int\s+main\s*\(", open(x,"r").read(), re.M) != None and
			re.search(r"//\s+do\s+not\s+build", open(x,"r").read(), re.M) == None)

gen_file_lists()

if use_verbose:
	if len(source_file_list)==0:
		print("No source files were found")
	else:
		print("Found the following source files:")
		print(" * "+"\n * ".join(source_file_list))
	print("")

	if len(header_file_list) == 0:
		print("No header files were found")
	else:
		print("Found the following header files:")
		print(" * "+"\n * ".join(header_file_list))
	print("")

	if len(main_file_list) == 0:
		print("No source file was identified as main file")
	else:
		print("The following source files were identified as main files:")
		print(" * "+"\n * ".join(main_file_list))
	print("")

################################################################################
## Generate Dependencies
##

def gen_dependencies():
	global all_header_include
	global local_header_include
	global direct_source_depend
	global link_against

	all_header_include = {}
	local_header_include = {}
	direct_source_depend = {}
	link_against = {}

	for x in source_file_list:
		all_include = find_all_includes(x)
		local_include = list(x for x in all_include if x in header_file_list)

		direct_depend = list(
			x for x in (
				os.path.splitext(f)[0] + ext
				for f in local_include
				for ext in source_extensions
			)
			if x in source_file_list
		)

		all_header_include[x] = all_include
		local_header_include[x] = local_include
		direct_source_depend[x] = direct_depend
		link_against[x] = direct_depend

	for x in source_file_list:
		visited = []
		unvisited = [x]
		while len(unvisited) != 0:
			y = unvisited.pop()
			visited.append(y)
			unvisited.extend(z for z in link_against[y] if not z in visited and not z in unvisited)
		link_against[x] = visited

gen_dependencies()

if use_verbose:
	for x in source_file_list:
		print("file "+x+" includes:")
		print(" * "+"\n * ".join(local_header_include[x]))
		print("directly depends on:")
		print(" * "+"\n * ".join(direct_source_depend[x]))
		print("and must be linked against:")
		print(" * "+"\n * ".join(link_against[x]))
		print("")

################################################################################
## Flag magic
##
def extract_flags(x):
	compile_self, compile_related, compile_all, link_all = extract_code_flags(x)
	h = list(os.path.split(y)[1] for y in all_header_include[x] if not y in header_file_list)

	if ("math.h" in h or "cmath" in h):
		link_all.append("-lm")

	if ("metis.h" in h):
		link_all.append("-lmetis")

	if ("kaffpa_interface.h" in h):
		link_all.append("-lkaffpa")

	if ("cl.h" in h or "opencl" in h):
		link_all.append("-lOpenCL")

	if ("omp.h" in h):
		link_all.append("-fopenmp")
		compile_self.append("-fopenmp")

	if ("thread" in h or "future" in h or "mutex" in h or "atomic" in h):
		link_all.append("-lpthread")

	return compile_self, compile_related, compile_all, link_all


def gen_flags():
	global compiler_flags
	global linker_flags

	compile_related = {}
	compile_all = {}
	link_all = {}

	compiler_flags = {}
	linker_flags = {}

	for x in source_file_list:
		compiler_flags[x] = []
		linker_flags[x] = []

	for x in source_file_list:
		compile_self, compile_related[x], compile_all[x], link_all[x] = extract_flags(x)
		compiler_flags[x].extend(compile_self)

	for x in source_file_list:
		compiler_flags[x].extend(compile_related[x])
		for y in direct_source_depend[x]:
			compiler_flags[x].extend(compile_related[y])
		for y in link_against[x]:
			compiler_flags[x].extend(compile_all[y])
			linker_flags[x].extend(link_all[y])

	for x in source_file_list:
		compiler_flags[x] = list(set(compiler_flags[x]))
		linker_flags[x] = list(set(linker_flags[x]))
gen_flags()

################################################################################
## Determine when source files were last modified
##


def gen_last_modified():
	global header_last_modify
	global source_last_modify
	header_last_modify = {}
	source_last_modify = {}
	for x in header_file_list:
		header_last_modify[x] = get_last_modified(x)
	for x in source_file_list:
		source_last_modify[x] = max(
			[get_last_modified(sys.argv[0]), get_last_modified(x)]+
			list(header_last_modify[y] for y in local_header_include[x]))

gen_last_modified()


################################################################################
## Build files
##

def build_all_files():

	for main_file in main_file_list:
		exe_file = os.path.splitext(main_file)[0]
		object_files = list(make_object_file_name(x) for x in link_against[main_file])

		if use_verbose:
			print("To build "+exe_file+" the following files must first be build:")
			print(" * "+"\n * ".join(object_files))
			print("")

		for x in link_against[main_file]:
			object_file = make_object_file_name(x)

			if get_last_modified(object_file) < source_last_modify[x]:
				if use_verbose:
					print("Compiling "+object_file)
				run_compiler(x, object_file, compiler_flags[x])
			else:
				if use_verbose:
					print("No need to recompile "+object_file)

		if get_last_modified(exe_file) < max([get_last_modified(x) for x in object_files]):
			if use_verbose:
				print("Linking "+exe_file)
			run_linker(object_files, exe_file, linker_flags[main_file])
		else:
			if use_verbose:
				print("No need to relink "+exe_file)

build_all_files()

################################################################################
## Cleanup
##
def do_cleanup():
	for x in (make_object_file_name(y) for y in source_file_list):
		if os.path.exists(x):
			cmd(["rm", x])

if use_clean:
	do_cleanup()

print("Build done")
