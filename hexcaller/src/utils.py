import datetime
import sys

def load_covs(filename):
    """Load coverage file
    args:
        filename: filename to load

    returns:
        dict: dictionary of samples, meta-data and raw coverage values
    """
    covs = {}
    linecount = 0

    with open(filename, 'r') as infile:
        for line in infile:
            dat = line.strip().split(',')
            linecount += 1
            if linecount <= 1:
                if dat[0] != "sample":
                    pass #log_error(f"invalid file format: {line} {dat[0]}")
            elif len(dat) == 7:
                if dat[0] in covs.keys():
                    covs[dat[0]].append({
                        'chr': dat[1],
                        'position': int(dat[2]),
                        'scoreA': float(dat[4]),
                        'scoreB': float(dat[6]),
                        'info': line.strip()
                    })
                else:
                    covs[dat[0]] = [{
                        'chr': dat[1],
                        'position': int(dat[2]),
                        'scoreA': float(dat[4]),
                        'scoreB': float(dat[6]),
                        'info': line.strip()
                    }]
    return covs


def log_error(message):
    """
    log error message with datetime and exit the program
    :param message:
    """

    dt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{dt}]: ERROR {message}")
    sys.exit(1)


def log_message(message):
    """
    log message
    :param message:
    """

    dt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sys.stderr.write(f"[{dt}]: {message}\n")
