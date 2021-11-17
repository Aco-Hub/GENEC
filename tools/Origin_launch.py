import sys
import os
import argparse
import shutil
import time
import socket
from six.moves import input
import platform
notify = False

if (
    platform.system() == 'Darwin'
    and int(platform.release()[:platform.release().find('.')]) > 10
):
    try:
        import pync
        notify = True
    except:
        print('pync not installed, no notifications on stops')


def mymail(
        email1,
        email2,
        message,
):
    with open('tmpf', 'w', encoding="utf8") as tmpfile:
        tmpfile.write(message)
    current_dir = os.getcwd()
    MyCommand = (
        'cat - tmpf << EOF | /usr/sbin/sendmail -r'
        + email1 + ' -t\nto:'
        + email2 + '\nsubject:'
        + current_dir + '\n\nEOF'
    )
    os.system(MyCommand)
    os.system('rm tmpf')


def stop_notif(
        current_dir,
        message,
):
    pync.notify(message, title=current_dir, sound='default')


# =======================================================================================
# User-specific default values are read in environment variables.
# If they do not exist, they are created and written in the file
# ~/.bash_profile.


def read_input_card(StarName):
    with open(
        f'{StarName}.input', 'r',
        encoding="utf8", errors="surrogateescape"
    ) as input_file:
        return input_file.read()

def write_input_card(StarName, input_card):
    with open(
        f'{StarName}.input', 'w',
        encoding="utf8", errors="surrogateescape"
    ) as input_file:
        return input_file.write(input_card)


def new_argument_parser(
        default_prog=None, email_address1=None, default_ncalc=1000,
):
    parser = argparse.ArgumentParser(
        description='Arguments for the launch of stellar models computation',
        usage=(
            f'''Origin_launch.py #star_name
--------------------------------------------------------
The default program is:
    {default_prog}
You can change it with option -e #prog

You can specify a phase (or model) to stop at with option -p #phase (-m #model)

More details on the options by calling Origin_launch.py -h
--------------------------------------------------------
'''
        )
    )
    parser.add_argument('StarName', help='Star name.', type=str)
    parser.add_argument(
        '-e', '--exe', help='Program to be used.', type=str,
        default=default_prog)
    parser.add_argument(
        '-p', '--phase', help='Phase to stop.', type=int)
    parser.add_argument(
        '-m', '--model', help='Model to stop.', type=int)
    parser.add_argument(
        '-i', '--initial', help='Initial file.', type=str)
    parser.add_argument(
        '-z', '--zip', action='store_true',
        help='Activate zipping all files after a series computation',
    )
    parser.add_argument(
        '-a', '--admail', help='mail address "To:"', type=str,
        default=email_address1
    )
    parser.add_argument(
        '-c', '--calcdir', help='Calculation directory.', type=str)
    parser.add_argument(
        '-n', '--ncalc', help='number of models computed before copying back',
        type=int, default=default_ncalc)
    parser.add_argument(
        '-N', '--NoMail', help='Activate Nomail mode', action='store_false')
    parser.add_argument(
        '-l', '--loop',
        help=(
            'Activate loop mode, up from top, down from bottom,'
            'upf or downf for ForceMode'
        ),
        type=str, default="")
    parser.add_argument(
        '-f', '--force',
        help='force looping, even in case of crash. Use with caution.',
        action='store_true')

    return parser.parse_args()


def main():

    # =======================================================================================
    source_dir = os.path.dirname(os.path.abspath(__file__))
    current_dir = os.getcwd()
    print('Current dir: '+current_dir)
    MakeInput = os.path.join(source_dir, 'MakeInput.py')
    # =======================================================================================
    loop_step = 0.
    loop_max = 0.033
    loop_min = 0.0005
    restart_loop = False
    time_to_transfer = False
    # =======================================================================================
    try:
        default_prog = os.environ['GENEC_DEFAULT_PROGRAM']
    except KeyError:
        default_prog = input(
            'You did not yet set the needed environment variables.\n'
            'Enter the default program to be used (full path): '
        )
        os.environ['GENEC_DEFAULT_PROGRAM'] = default_prog

        with open(
            os.path.expanduser('~/.bash_profile'), 'a',
            encoding="utf8", errors="surrogateescape",
        ) as bp:
            bp.write(
                '##################\n'
                '# genec variables \n'
                '##################\n'
            )
            bp.write('export GENEC_DEFAULT_PROGRAM="'+default_prog+'"\n')
    try:
        email_address1 = os.environ['GENEC_EMAIL_ADDRESS']
        print('email address for sender: '+email_address1)
    except KeyError:
        email_address1 = input('Enter the email address to be used (sender): ')
        os.environ['GENEC_EMAIL_ADDRESS'] = email_address1
        with open(
            os.path.expanduser('~/.bash_profile'), 'a',
            encoding="utf8", errors="surrogateescape",
        ) as bp:
            bp.write(f'export GENEC_EMAIL_ADDRESS="{email_address1}"'+'\n')

    # Defaults
    default_ncalc = 1000

    args = new_argument_parser(
        default_prog=default_prog,
        default_ncalc=default_ncalc,
        email_address1=email_address1,
    )

    ProgEvol = args.exe
    StarName = args.StarName
    phase_stop = args.phase
    model_stop = args.model
    initial_file = args.initial
    email_address2 = args.admail
    calc_dir = args.calcdir
    ncalc = args.ncalc
    MailMode = args.NoMail
    Zipping = args.zip
    LoopMode = args.loop
    ForceMode = args.force
    if LoopMode and LoopMode[-1] == 'f':
        ForceMode = True
    # =======================================================================================
    if platform.system() == 'Darwin':
        MailMode = False
    base_mailmode = MailMode

    if "down" in LoopMode:
        initial_loop = [0.0005, 0.005]
        loop_step = 0.0005
    elif "up" in LoopMode:
        initial_loop = [0.033, 0.033]
        loop_step = -0.0005
    print(loop_step)

    if calc_dir is None:
        calc_dir = ''
    if initial_file is None:
        initial_file = ''

    modanfs = 'modanf='
    nwseqs = 'nwseq='
    nzmods = 'nzmod='
    phases = 'phase='
    xcns = 'xcn='
    deltal = 'deltal='
    deltat = 'deltat='
    CommandLaunch = ProgEvol+' < '+StarName+'.input'

    print('Prog: '+ProgEvol)
    print('StarName: '+StarName)
    if calc_dir:
        print('Calc dir: '+calc_dir)
    if initial_file:
        print('starting on initial file: '+initial_file)
    if calc_dir != '':
        Zipping = True
    if initial_file == '':
        Inputs = read_input_card(StarName)
        ibfile = Inputs.rfind(modanfs)+len(modanfs)
        ibfile_end = Inputs[ibfile:].find('\n')
        imod = Inputs.rfind(nwseqs)+len(nwseqs)
        imod_end = Inputs[imod:].find('\n')
        iphase = Inputs.rfind(phases)+len(phases)
        iphase_end = Inputs[iphase:].find('\n')
        modanf = int(Inputs[ibfile:ibfile+ibfile_end])
        nwseq = int(Inputs[imod:imod+imod_end])
        phase = int(Inputs[iphase:iphase+iphase_end])

    mytime = time.strftime('%A %d %B %Y at %H:%M:%S')
    with open(
        'computation.log', 'a',
        encoding="utf8", errors="surrogateescape",
    ) as logfile:
        logfile.write(
            f'Computation on {socket.gethostname()} started on {mytime}'+'\n')
        logfile.write(f'Directory: {current_dir}'+'\n')
        if initial_file == '':
            logfile.write(
                f'starting model is {nwseq}, phase: {phase}'+'\n')
        if calc_dir != '':
            logfile.write(
                f'Computation performed in directory: {calc_dir}'+'\n')
        logfile.write(f'Program used: {ProgEvol}'+'\n')
        if initial_file != '':
            logfile.write(
                f'Program launched on initial file {initial_file}'+'\n')
        if phase_stop is not None:
            logfile.write(
                f'Requested stop at phase {phase_stop}'+'\n')
        if model_stop is not None:
            logfile.write(
                'Requested stop at model {model_stop}'+'\n')

    time_start = time.time()
    answer = ''

    answer = ''

    if initial_file != '':
        if StarName+'.input' in os.listdir('.'):
            answer = input(
                'This star seems to be already partially computed.\n'
                f'Are you sure you want to proceed from file {initial_file}?'
                '\n yes(y) or no(n): '
            )
            if not answer:
                sys.exit()
        if answer.lower() in 'yes':
            os.system(ProgEvol+' < '+initial_file)
            try:
                runlog = open(
                    'runfile', 'r', encoding="utf8", errors="surrogateescape"
                )
            except IOError:
                sys.exit()
            runstat = runlog.read().strip(' \n\t')
            if runstat != 'running':
                if runstat != '':
                    print('Program stopped with message: '+runstat)
                else:
                    print('Program aborted...')
                sys.exit()
            else:
                if Zipping:
                    CommandZip = (
                        f'gzip -f {StarName}.[l,v,x,y]0000001'
                        f' {StarName}.b00000 {StarName}.b00001'
                        f' {StarName}_StrucData_0000001.dat'
                    )
                    os.system(CommandZip)

    relaunch_advection = [True, 0, 0]

    if calc_dir != '':
        try:
            os.makedirs(calc_dir)
            print('Calculation directory '+calc_dir+' successfully created')
        except OSError:
            print('Calculation directory '+calc_dir+' already exists')
        try:
            os.remove(os.path.join(calc_dir, '*'))
        except:
            pass

    # Main loop
    while True:
        Inputs = read_input_card(StarName)
        ibfile = Inputs.rfind(modanfs)+len(modanfs)
        ibfile_end = Inputs[ibfile:].find('\n')
        imod = Inputs.rfind(nwseqs)+len(nwseqs)
        imod_end = Inputs[imod:].find('\n')
        iphase = Inputs.rfind(phases)+len(phases)
        iphase_end = Inputs[iphase:].find('\n')
        modanf = int(Inputs[ibfile:ibfile+ibfile_end])
        nwseq = int(Inputs[imod:imod+imod_end])
        phase = int(Inputs[iphase:iphase+iphase_end])
        if calc_dir != '':
            if nwseq % ncalc == 1:
                time_to_transfer = True
            else:
                time_to_transfer = False
        if LoopMode != "" and restart_loop:
            MailMode = False
            LineLeft = Inputs.rfind(deltal)
            LineRight = Inputs.rfind(deltat) + 14
            NewLine = (
                f'{deltal}{initial_loop[0]:0>7.5f}, '
                f'{deltat}{initial_loop[1]:0>7.5f}'
            )
            ToBeReplaced = Inputs[LineLeft:LineRight]
            Inputs = Inputs.replace(ToBeReplaced, NewLine)
            with open(
                StarName+'.input', 'w',
                encoding="utf8", errors="surrogateescape",
            ) as InputFile:
                InputFile.write(Inputs)
            if "up" in LoopMode:
                if initial_loop[1] <= loop_min:
                    initial_loop[0] = initial_loop[0]+loop_step
                    initial_loop[1] = loop_max
                else:
                    initial_loop[1] = initial_loop[1]+loop_step
            if "down" in LoopMode:
                if initial_loop[1] >= loop_max:
                    initial_loop[0] = initial_loop[0]+loop_step
                    initial_loop[1] = loop_min
                else:
                    initial_loop[1] = initial_loop[1]+loop_step

        if phase_stop is not None and phase == phase_stop:
            stop_message = f'{nwseq} : Phase {phase_stop} reached.'
            print(stop_message)
            if MailMode and len(email_address2) != 0:
                mymail(email_address1, email_address2, stop_message)
            elif notify is True:
                stop_notif(current_dir, stop_message)
            break

        if model_stop is not None and nwseq > model_stop:
            stop_message = f'Model {model_stop} reached.'
            print(stop_message)
            if MailMode and len(email_address2) != 0:
                mymail(email_address1, email_address2, stop_message)
            elif notify is True:
                stop_notif(current_dir, stop_message)
            break

        print(f'New model {nwseq} with bfile {modanf}')
        print(f'Phase: {phase}')

        needed_net = [i for i in os.listdir('.') if i[0:3] == 'net']
        needed_for_calc = needed_net+['input_changes.log', f'{StarName}.input']
        if os.path.isfile(
            f'{StarName}.b{modanf:05d}.gz'
        ):
            needed_for_calc = (
                needed_for_calc + [f'{StarName}.b{modanf:05d}.gz']
            )
        else:
            needed_for_calc = (
                needed_for_calc + [f'{StarName}.b{modanf:05d}']
            )
        try:
            needed_for_calc = needed_for_calc+[f'.PlotData_{StarName}']
        except:
            pass

        if calc_dir != '' and os.getcwd() != calc_dir:
            for file in needed_for_calc:
                shutil.copy2(file, calc_dir)
            os.chdir(calc_dir)

        try:
            os.remove('runfile')
        except:
            pass
        if Zipping:
            if os.path.isfile(f'{StarName}.b{modanf:05d}.gz'):
                os.system(f'gunzip {StarName}.b{modanf:05d}.gz')

        # actual call to GENEC
        os.system(CommandLaunch)

        if Zipping:
            CommandZip = (
                f'gzip -f {StarName}.[l,v,x,y]{nwseq:07d} '
                f'{StarName}.b{modanf:05d} {StarName}.b{modanf+1:05d} '
                f'{StarName}_StrucData_{nwseq:07d}.dat'
            )
            os.system(CommandZip)

        if calc_dir != '' and time_to_transfer:
            print(f'time to transfer ({nwseq})')
            result_files = [
                i for i in os.listdir('.') if i[0:4] == StarName[0:4]
            ]
            result_files = result_files+['input_changes.log']
            try:
                result_files = result_files+['.PlotData_'+StarName]
            except:
                pass
            for file in result_files:
                shutil.move(file, current_dir+'/'+file)
            shutil.copy2('runfile', current_dir)
            os.chdir(current_dir)

        try:
            runlog = open('runfile', 'r', encoding="utf8", errors="surrogateescape")
        except:
            if (
                "up" in LoopMode
                and (
                    initial_loop[0] >= loop_min
                    or initial_loop[1] >= loop_min
                )
            ):
                restart_loop = True
            elif (
                "down" in LoopMode
                and (
                    initial_loop[0] <= loop_max
                    or initial_loop[1] <= loop_max
                )
            ):
                restart_loop = True
            else:
                stop_message = f'Program stopped with message: {runstat}'
                if MailMode and len(email_address2) != 0:
                    mymail(email_address1, email_address2, runstat)
                elif notify is True:
                    stop_notif(current_dir, stop_message)
                break
            continue
        runstat = runlog.read().strip(' \n\t')
        if runstat != 'running':
            if runstat != '':
                if 'ZAMS reached' in runstat:
                    input_card = read_input_card(StarName)
                    nwseq = int(
                        input_card[
                            input_card.rfind(nwseqs)
                            + len(nwseqs):input_card.find('\n modanf')
                        ]
                    )
                    modanf = int(
                        input_card[
                            input_card.rfind(modanfs)
                            + len(modanfs):input_card.find('\n nzmod')
                        ]
                    )
                    nzmod = int(
                        input_card[
                            input_card.rfind(nzmods)
                            + len(nzmods):input_card.find('\n&END')
                        ]
                    )
                    print(
                        f'ZAMS reached: NWSEQ= {nwseq}, MODANF= {modanf}, '
                        f'NZMOD= {nzmod}'
                    )
                    if Zipping:
                        if os.path.isfile(f'{StarName}.b{modanf:05d}.gz'):
                            os.system(f'gunzip {StarName}.b{modanf:05d}.gz')
                    os.system(CommandLaunch)
                    if Zipping:
                        CommandZip = (
                            f'gzip -f {StarName}.[l,v,x,y]{nwseq:07d} '
                            f'{StarName}.b{modanf:05d} '
                            f'{StarName}_StrucData_{nwseq:07d}.dat'
                        )
                        os.system(CommandZip)
                    with open('runfile', 'r') as newrun:
                        runlog = newrun.read().strip(' \n\t')
                        if runlog == 'running':
                            input_card = read_input_card(StarName)
                            nzmod = (
                                int(
                                    input_card[
                                        input_card.rfind(nzmods)
                                        + len(nzmods):input_card.find('\n&END')
                                    ]
                                )
                            )
                            input_card = input_card.replace(
                                'nzmod='+str(nzmod), 'nzmod=10')
                            input_card = input_card.replace(
                                'gkorm=.300', 'gkorm=.100')
                            with open(
                                StarName+'.input', 'w', encoding="utf8",
                                errors="surrogateescape",
                            ) as input_file:
                                input_file.write(input_card)
                        else:
                            break
                elif (
                    'Problem during advection' in runstat
                    or 'Advection not applied' in runstat
                    or 'Problem with conservation of angular momentum during advection' in runstat
                    or 'Ang. mom. variation too large during diffusion' in runstat
                ) and (LoopMode[-1] != 'f' or not ForceMode):
                    timestep = int(runstat[0:runstat.find(':')])
                    if relaunch_advection[2] == timestep:
                        stop_message = 'Program stopped with message: '+runstat
                        print('Automatic relaunch for advection failed at previous sequence already twice. ')
                        print(' Please retry with a smaller time step.')
                        if MailMode and len(email_address2) != 0:
                            mymail(email_address1, email_address2, runstat)
                        elif notify is True:
                            stop_notif(current_dir, runstat)
                        break
                    relaunch_advection[2] = timestep
                    if relaunch_advection[0]:
                        input_card = read_input_card(StarName)
                        nwseq = int(
                            input_card[
                                input_card.rfind(nwseqs)
                                + len(nwseqs):input_card.find('\n modanf')
                            ]
                        )
                        if timestep % 10 == 1:
                            MakeCommand = MakeInput + " " + StarName + " " + str(nwseq-10)
                            os.system(MakeCommand)
                        input_card = read_input_card(StarName)
                        to_replace = (
                            "xcn="
                            + input_card[
                                input_card.rfind(xcns)
                                + len(xcns):input_card.find(
                                    '\n&END', input_card.rfind(xcns)+len(xcns))
                            ]
                        )
                        xcn = float(
                            input_card[
                                input_card.rfind(xcns)
                                + len(xcns):input_card.find(
                                    '\n&END', input_card.rfind(xcns)+len(xcns))
                            ]
                        )
                        input_card = input_card.replace(
                            to_replace, 'xcn=0.300')
                        write_input_card(StarName, input_card)
                        relaunch_advection[1] = relaunch_advection[1] + 1
                        if timestep % 10 == 1 or relaunch_advection[1] > 1:
                            relaunch_advection[0] = False
                    else:
                        stop_message = 'Program stopped with message: '+runstat
                        if MailMode and len(email_address2) != 0:
                            mymail(email_address1, email_address2, runstat)
                        elif notify is True:
                            stop_notif(current_dir, runstat)
                        break
                else:
                    if (
                        "up" in LoopMode
                            and (
                                initial_loop[0] >= loop_min
                                or initial_loop[1] >= loop_min
                            )
                    ):
                        restart_loop = True
                    elif (
                        "down" in LoopMode
                        and (
                            initial_loop[0] <= loop_max
                            or initial_loop[1] <= loop_max
                        )
                    ):
                        restart_loop = True
                    else:
                        stop_message = 'Program stopped with message: '+runstat
                        if MailMode and len(email_address2) != 0:
                            mymail(email_address1, email_address2, runstat)
                        elif notify is True:
                            stop_notif(current_dir, runstat)
                        break
            elif ForceMode:
                if (
                    "up" in LoopMode
                        and (
                            initial_loop[0] >= loop_min
                            or initial_loop[1] >= loop_min
                        )
                ):
                    restart_loop = True
                elif (
                    "down" in LoopMode
                    and (
                        initial_loop[0] <= loop_max
                        or initial_loop[1] <= loop_max
                    )
                ):
                    restart_loop = True
                else:
                    stop_message = 'Program stopped with message: '+runstat
                    if MailMode and len(email_address2) != 0:
                        mymail(email_address1, email_address2, runstat)
                    elif notify is True:
                        stop_notif(current_dir, runstat)
                    break
                continue
            else:
                stop_message = 'Program aborted...'
                break
        else:
            restart_loop = False
            MailMode = base_mailmode
            if "up" in LoopMode:
                initial_loop = [loop_max, loop_max]
            else:
                initial_loop = [loop_min, loop_min]

        if (
            'Problem during advection' not in runstat
            and 'Advection not applied' not in runstat
            and 'Problem with conservation of angular momentum during advection' not in runstat
        ):
            relaunch_advection[0] = True
            relaunch_advection[1] = 0

    if calc_dir != '' and os.getcwd() != current_dir:
        result_files = [i for i in os.listdir('.') if i[0:4] == StarName[0:4]]
        result_files = result_files+['input_changes.log']
        try:
            result_files = result_files+['.PlotData_'+StarName]
        except:
            pass
        for file in result_files:
            shutil.move(file, current_dir+'/'+file)
        shutil.copy2('runfile', current_dir)
        os.chdir(current_dir)
    time_stop = time.time()
    diff_time = int(round(time_stop-time_start))
    m, s = divmod(diff_time, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    computation_time = (
        f'{d:d} day(s) {h:d} hour(s) {m:d} minutes and {s:d} seconds'
    )
    with open(
        'computation.log', 'a',
        encoding="utf8",  errors="surrogateescape"
    ) as logfile:
        logfile.write(stop_message+'\n')
        logfile.write('Computation lasted '+computation_time+'\n')
        logfile.write('----------------------------------------\n\n')


if __name__ == "__main__":
    main()
