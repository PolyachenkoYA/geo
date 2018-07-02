################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/comp_geo.cu 

OBJS += \
./src/comp_geo.o 

CU_DEPS += \
./src/comp_geo.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-9.2/bin/nvcc -include ../../geo.cu -G -g -O0 -Xcompiler -fopenmp -std=c++11   -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-9.2/bin/nvcc -include ../../geo.cu -G -g -O0 -Xcompiler -fopenmp -std=c++11 --compile --relocatable-device-code=false  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


