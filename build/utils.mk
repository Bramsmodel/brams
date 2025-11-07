# gcc info
export GCC_VER_TO_INT = $(shell echo $(1) | awk -F. '{print $$3+100*($$2+100*$$1)}')
export GCC_VERSION := $(shell gcc -dumpfullversion)
export GCC_MAJOR_VERSION := $(firstword $(subst ., , $(GCC_VERSION)))
export GCC_VERSION_INT := $(call GCC_VER_TO_INT, $(GCC_VERSION))
export GCC_VER_GTE10 := $(shell [ $(GCC_VERSION_INT) -ge 100000 ] || echo 0 )

