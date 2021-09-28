void exponential_generator_int64_ (mypair<uint64_t, uint64_t>* A,
                                  int n, int exp_cutoff, double exp_lambda) {
    pbbs::sequence<uint32_t> nums(exp_cutoff);
    pbbs::sequence<mypair<uint64_t, uint64_t>> B(n);
    // double base = 1 - exp(-1);
    // cout << "1 - e^-1 = " << base << endl;
    // cout << "e^-2 = " << exp(-2) << endl;

    /* 1. making nums[] array */
    parallel_for (0, exp_cutoff, [&] (int i) {
    // for (int i = 0; i < exp_cutoff; i++) {
        nums[i] = (double)n * (exp(-exp_lambda * i) * (1 - exp(-exp_lambda)));
        // cout << nums[i] << endl;
    // }
    }, 1);

    double offset = pbbs::reduce(nums, pbbs::addm<uint32_t>());
    nums[0] += (n - offset);
    // cout << "offset/n = " << offset << "/" << n << endl;
    // checking if the sum of nums[] equals to n
    if (pbbs::reduce(nums, pbbs::addm<uint32_t>()) == (uint32_t)n) {
        cout << "sum of nums[] == n" << endl;
    }

    /* 2. scan to calculate position */
    uint32_t* addr = new uint32_t[exp_cutoff];
    parallel_for (0, exp_cutoff, [&] (uint32_t i) {
        addr[i] = nums[i];
    }, 1);
    scan_inplace__(addr, exp_cutoff); // store all addresses into addr[]

    /* 3. distribute random numbers into A[i].first */
    parallel_for (0, exp_cutoff, [&] (size_t i) {
        size_t st = (i == 0) ? 0 : addr[i-1],
               ed = (i == (uint32_t)exp_cutoff-1) ? n : addr[i];
        for (size_t j = st; j < ed; j++) {
            B[j].first = pbbs::hash64_2(i);
        }
    }, 1);
    parallel_for (0, n, [&] (size_t i){
        B[i].second = pbbs::hash64_2(i);
    }, 1);

    /* 4. shuffle the keys */
    pbbs::sequence<mypair<uint64_t, uint64_t>> C = pbbs::random_shuffle(B, n);

    parallel_for (0, n, [&] (size_t i) {
        A[i] = C[i];
    });

    delete[] addr;
}

// distributions
    if (strcmp(argv[2], "uniform") == 0) {
        uint64_t uniform_max_range = atoi(argv[3]);
        cout << "Uniform distribution (64bit key, 64bit value)..." << endl;
        cout << "Uniform parameter = " << uniform_max_range << endl;
        uniform_generator_int64_(A, n, uniform_max_range);
    } else if (strcmp(argv[2], "exponential") == 0) {
        double exp_lambda = stod(argv[3]);
        cout << "Exponential distribution (64bit key, 64bit value)..." << endl;
        cout << "Exponential parameter = " << exp_lambda << endl;
        exponential_generator_int64_(A, n, 1000000, exp_lambda);
    } else if (strcmp(argv[2], "zipfian") == 0) {
        uint32_t zipfian_s = stod(argv[3]);
        cout << "Zipfian distribution (64bit key, 64bit value)..." << endl;
        cout << "Zipfian parameter = " << zipfian_s << endl;
        zipfian_generator_int64_(A, n, zipfian_s);
    }