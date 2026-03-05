from concurrent.futures import ThreadPoolExecutor


def run_parallel(function, items, threads=4):

    results = []

    with ThreadPoolExecutor(max_workers=threads) as executor:

        futures = [executor.submit(function, item) for item in items]

        for f in futures:
            results.append(f.result())

    return results
