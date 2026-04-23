FC = gfortran
FFLAGS = -Wall -fimplicit-none -J$(BUILDDIR) -cpp -fallow-argument-mismatch
SRCDIR = ../src
BUILDDIR = ../build
BINDIR = ../bin

# Source modules from GENEC needed for tests (order matters)
GENEC_SRCS = \
	$(SRCDIR)/io_definitions.f90 \
	$(SRCDIR)/evol.f90 \
	$(SRCDIR)/storage.f90 \
	$(SRCDIR)/safestop.f90 \
	$(SRCDIR)/const.f90 \
	$(SRCDIR)/caramodele.f90 \
	$(SRCDIR)/inputparam.f90 \
	$(SRCDIR)/interpolation.f90 \
	$(SRCDIR)/nagmod.f90 \
	$(SRCDIR)/SmallFunc.f90 \
	$(SRCDIR)/strucmod.f90 \
	$(SRCDIR)/ionisation.f90

GENEC_OBJS = \
	$(BUILDDIR)/io_definitions.o \
	$(BUILDDIR)/evol.o \
	$(BUILDDIR)/storage.o \
	$(BUILDDIR)/safestop.o \
	$(BUILDDIR)/const.o \
	$(BUILDDIR)/caramodele.o \
	$(BUILDDIR)/inputparam.o \
	$(BUILDDIR)/interpolation.o \
	$(BUILDDIR)/nagmod.o \
	$(BUILDDIR)/SmallFunc.o \
	$(BUILDDIR)/strucmod.o \
	$(BUILDDIR)/ionisation.o

# Test framework
FRAMEWORK_OBJ = $(BUILDDIR)/test_framework.o

# Test executables
TESTS = \
	$(BINDIR)/test_interpolation \
	$(BINDIR)/test_smallfunc \
	$(BINDIR)/test_const \
	$(BINDIR)/test_ionisation

.PHONY: all clean run bench

all: dirs $(TESTS)

bench: dirs $(GENEC_OBJS) $(BINDIR)/bench_fortran

dirs:
	mkdir -p $(BUILDDIR) $(BINDIR)

# Build GENEC modules
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Build test framework
$(BUILDDIR)/test_framework.o: test_framework.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Build and link test executables
$(BINDIR)/test_interpolation: test_interpolation.f90 $(FRAMEWORK_OBJ) $(GENEC_OBJS)
	$(FC) $(FFLAGS) -o $@ $(BUILDDIR)/evol.o $(BUILDDIR)/interpolation.o $(FRAMEWORK_OBJ) $<

$(BINDIR)/test_smallfunc: test_smallfunc.f90 $(FRAMEWORK_OBJ) $(GENEC_OBJS)
	$(FC) $(FFLAGS) -o $@ $(BUILDDIR)/io_definitions.o $(BUILDDIR)/evol.o $(BUILDDIR)/storage.o $(BUILDDIR)/safestop.o $(BUILDDIR)/const.o $(BUILDDIR)/caramodele.o $(BUILDDIR)/inputparam.o $(BUILDDIR)/interpolation.o $(BUILDDIR)/nagmod.o $(BUILDDIR)/SmallFunc.o $(FRAMEWORK_OBJ) $<

$(BINDIR)/test_const: test_const.f90 $(FRAMEWORK_OBJ) $(GENEC_OBJS)
	$(FC) $(FFLAGS) -o $@ $(BUILDDIR)/evol.o $(BUILDDIR)/const.o $(FRAMEWORK_OBJ) $<

$(BINDIR)/test_ionisation: test_ionisation.f90 $(FRAMEWORK_OBJ) $(GENEC_OBJS)
	$(FC) $(FFLAGS) -o $@ $(BUILDDIR)/io_definitions.o $(BUILDDIR)/evol.o $(BUILDDIR)/storage.o $(BUILDDIR)/safestop.o $(BUILDDIR)/const.o $(BUILDDIR)/strucmod.o $(BUILDDIR)/ionisation.o $(FRAMEWORK_OBJ) $<

run: all
	@echo ""
	@echo "Running all tests..."
	@echo "==============================="
	@echo "Test,Time_ms" > test_benchmarks.csv
	@failed=0; \
	total_start=$$(date +%s%N); \
	for test in $(TESTS); do \
		echo ""; \
		echo ">>> Running $$(basename $$test)..."; \
		start=$$(date +%s%N); \
		$$test || failed=$$((failed + 1)); \
		end=$$(date +%s%N); \
		elapsed=$$(( (end - start) / 1000000 )); \
		echo "    Time: $${elapsed}ms"; \
		echo "$$(basename $$test),$${elapsed}" >> test_benchmarks.csv; \
	done; \
	total_end=$$(date +%s%N); \
	total=$$(( (total_end - total_start) / 1000000 )); \
	echo ""; \
	echo "==============================="; \
	echo "Total time: $${total}ms"; \
	echo "Total,$${total}" >> test_benchmarks.csv; \
	echo "Results saved to test_benchmarks.csv"; \
	if [ $$failed -eq 0 ]; then \
		echo "All test suites passed."; \
	else \
		echo "$$failed test suite(s) had failures."; \
		exit 1; \
	fi

# Build benchmark harness
$(BINDIR)/bench_fortran: bench_fortran.f90 $(GENEC_OBJS)
	$(FC) $(FFLAGS) -O3 -o $@ $(BUILDDIR)/io_definitions.o $(BUILDDIR)/evol.o $(BUILDDIR)/storage.o $(BUILDDIR)/safestop.o $(BUILDDIR)/const.o $(BUILDDIR)/caramodele.o $(BUILDDIR)/inputparam.o $(BUILDDIR)/interpolation.o $(BUILDDIR)/nagmod.o $(BUILDDIR)/SmallFunc.o $(BUILDDIR)/strucmod.o $(BUILDDIR)/ionisation.o $<

clean:
	rm -f $(BUILDDIR)/test_*.o
	rm -f $(TESTS) $(BINDIR)/bench_fortran
