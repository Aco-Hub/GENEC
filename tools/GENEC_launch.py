"""
Launcher for GENEC
"""
import sys
import os
import argparse
import shutil
import time
import socket
import platform
from six.moves import input

PROGRAM_NAME = os.path.basename(__file__)
ENCODING = 'utf8'
NOTIFY = False

if (
    platform.system() == 'Darwin'
    and int(platform.release()[:platform.release().find('.')]) > 10
):
    try:
        import pync
        NOTIFY = True
    except ModuleNotFoundError:
        print('pync not installed, no notifications on stops')


def mymail(email1, email2, message):
    "send email notification"
    with open('tmpf', 'w', encoding=ENCODING) as tmpfile:
        tmpfile.write(message)
    current_dir = os.getcwd()
    my_command = (
        f'cat - tmpf << EOF | /usr/sbin/sendmail -r{email1} -t\n'
        f'to:{email2}\n'
        f'subject:{current_dir}\n'
        f'\n'
        f'EOF'
    )
    os.system(my_command)
    os.system('rm tmpf')


def stop_notify(current_dir, message):
    "send notification via pync"
    pync.notify(message, title=current_dir, sound='default')


def new_argument_parser(genec_defaults):
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Arguments for the launch of stellar models computation',
        usage=(
            f'{PROGRAM_NAME} #star_name\n'
            f'--------------------------------------------------------\n'
            f'The default program is:\n'
            f'   {genec_defaults["program"]}\n'
            f'You can change it with option -e #prog\n'
            f'\n'
            f'You can specify a phase (or model) to stop at with option '
            f'-p #phase (-m #model)\n'
            f'\n'
            f'More details on the options by calling {PROGRAM_NAME} -h\n'
            f'--------------------------------------------------------\n'
        )
    )

    parser.add_argument(
        'StarName', dest='star_name',
        help='Star name', type=str
    )
    parser.add_argument(
        '-e', '--exe',
        dest='program',
        help='Program to be used.', type=str,
        default=genec_defaults['program']
    )
    parser.add_argument(
        '-p', '--phase', dest='phase_stop',
        help='Phase to stop.', type=int
    )
    parser.add_argument(
        '-m', '--model', dest='model_stop',
        help='Model to stop.', type=int
    )
    parser.add_argument(
        '-i', '--initial', dest='initial_file',
        help='Initial file.', type=str
    )
    parser.add_argument(
        '-z', '--zip', dest='zipping',
        help='Activate zipping all files after a series computation',
        action='store_true'
    )
    parser.add_argument(
        '-a', '--admail', dest='email_receiver',
        help='mail address "To:"', type=str,
        default=genec_defaults['email_sender']
    )
    parser.add_argument(
        '-c', '--calcdir', dest='calc_dir',
        help='Calculation directory.', type=str
    )
    parser.add_argument(
        '-n', '--ncalc',
        help='number of models computed before copying back', type=int,
        default=genec_defaults['ncalc']
    )
    parser.add_argument(
        '-N', '--NoMail', dest='mail_mode',
        help='Activate no e-mail mode', action='store_false'
    )
    parser.add_argument(
        '-l', '--loop', dest='loop_mode',
        help=(
            'Activate loop mode, up from top, down from bottom, '
            'upf or downf for force_mode'
        ),
        type=str, default=""
    )
    parser.add_argument(
        '-f', '--force', dest='force_mode',
        help='force looping, even in case of crash. Use with caution.',
        action='store_true'
    )

    return parser.parse_args()


def main():
    """
    Main GENEC runner
    """
    # =======================================================================================
    # User-specific default values are read in environment variables.
    # If they do not exist, they are created and written in the file
    # ~/.bash_profile.
    genec_defaults = {}
    try:
        genec_defaults['program'] = os.environ['GENEC_DEFAULT_PROGRAM']
    except KeyError:
        if os.path.isfile('../bin/genec'):
            use_default = input(
                'You did not yet set the needed environment variables.\n'
                'Using the default "../bin/genec"\n'
                'Is this ok (Y/N)? '
            )
            if use_default.lower() == "y":
                genec_defaults['program'] = "../bin/genec"
            else:
                genec_defaults['program'] = input(
                    'Enter the program to be used (full path): '
                )
        write_default = input(
            'Write the default program to your bash profile (Y/N)? '
        )
        if write_default.lower() == "y":
            os.environ['GENEC_DEFAULT_PROGRAM'] = genec_defaults['program']
            with open(
                os.path.expanduser('~/.bash_profile'), 'a', encoding=ENCODING,
            ) as bp:
                bp.write(
                    '##################\n'
                    '# genec variables \n'
                    '##################\n'
                )
                bp.write(
                    f'export GENEC_DEFAULT_PROGRAM="'
                    f'{genec_defaults["program"]}"\n'
                )
    try:
        genec_defaults['email_sender'] = os.environ['GENEC_EMAIL_ADDRESS']
        print(f'email address for sender: {genec_defaults["email_sender"]}')
    except KeyError:
        genec_defaults['email_sender'] = input(
            'Enter the email address to be used (sender): '
        )
        os.environ['GENEC_EMAIL_ADDRESS'] = genec_defaults['email_sender']
        with open(
                os.path.expanduser('~/.bash_profile'), 'a', encoding=ENCODING
        ) as bp:
            bp.write(
                f'export GENEC_EMAIL_ADDRESS="'
                f'{genec_defaults["email_sender"]}"\n'
            )

    # =======================================================================================
    source_dir = os.path.dirname(os.path.abspath(__file__))
    current_dir = os.getcwd()
    print(f'Current dir: {current_dir}')
    make_input = os.path.join(source_dir, 'MakeInput.py')
    # =======================================================================================
    loop_step = 0.
    loop_max = 0.033
    loop_min = 0.0005
    restart_loop = False
    genec_defaults['ncalc'] = 1000
    time_to_transfer = False
    requested_stop = False
    # =======================================================================================

    settings = vars(new_argument_parser(genec_defaults))

    program = settings['program']
    star_name = settings['star_name']
    phase_stop = settings['phase_stop']
    model_stop = settings['model_stop']
    initial_file = settings['initial_file']
    email_receiver = settings['email_receiver']
    calc_dir = settings['calc_dir']
    ncalc = settings['ncalc']
    mail_mode = settings['mail_mode']
    zipping = settings['zipping']
    loop_mode = settings['loop_mode']
    force_mode = settings['force_mode']
    if loop_mode and loop_mode[-1] == 'f':
        force_mode = True
    # =======================================================================================
    if platform.system() == 'Darwin':
        mail_mode = False
    base_mail_mode = mail_mode

    if "down" in loop_mode:
        initial_loop = [0.0005, 0.005]
        loop_step = 0.0005
    elif "up" in loop_mode:
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
    endphases = 'end_at_phase='
    endmodels = 'end_at_model='
    phases = 'phase='
    xcns = 'xcn='
    deltal = 'deltal='
    deltat = 'deltat='
    command_launch = f'{program} < {star_name}.input'

    print(f'Prog: {program}')
    print(f'StarName: {star_name}')
    if calc_dir:
        print(f'Calc dir: {calc_dir}')
    if initial_file:
        print('starting on initial file: {initial_file}')
    if calc_dir != '':
        zipping = True
    if initial_file == '':
        with open(f'{star_name}.input', 'r', encoding=ENCODING) as input_file:
            inputs = input_file.read()
            ibfile = inputs.rfind(modanfs)+len(modanfs)
            ibfile_end = inputs[ibfile:].find('\n')
            imod = inputs.rfind(nwseqs)+len(nwseqs)
            imod_end = inputs[imod:].find('\n')
            iphase = inputs.rfind(phases)+len(phases)
            iphase_end = inputs[iphase:].find('\n')
            modanf = int(inputs[ibfile:ibfile+ibfile_end])
            nwseq = int(inputs[imod:imod+imod_end])
            phase = int(inputs[iphase:iphase+iphase_end])

    mytime = time.strftime('%A %d %B %Y at %H:%M:%S')
    with open('computation.log', 'a', encoding=ENCODING) as logfile:
        logfile.write(
            f'Computation on {socket.gethostname()} started on {mytime}\n'
        )
        logfile.write('Directory: {current_dir}\n')
        if initial_file == '':
            logfile.write(
                f'starting model is {str(nwseq)}, phase: {str(phase)}\n'
            )
        if calc_dir != '':
            logfile.write(f'Computation performed in directory: {calc_dir}\n')
        logfile.write(f'Program used: {program}\n')
        if initial_file != '':
            logfile.write(f'Program launched on initial file {initial_file}\n')
        if phase_stop is not None:
            logfile.write(f'Requested stop at phase {str(phase_stop)}\n')
        if model_stop is not None:
            logfile.write(f'Requested stop at model {str(model_stop)}\n')

    time_start = time.time()

    answer = ''

    if initial_file != '':
        if f'{star_name}.input' in os.listdir('.'):
            answer = input(
                f'This star seems to be already partially computed.\n'
                f'Are you sure you want to proceed from file {initial_file}?\n'
                f' yes(y) or no(n): '
            )
            if not answer:
                sys.exit()
        if answer.lower() in 'yes':
            os.system(f'{program} < {initial_file}')
            try:
                runlog = open('runfile', 'r', encoding=ENCODING)
            except OSError:
                sys.exit()
            runstat = runlog.read().strip(' \n\t')
            if runstat != 'running':
                if runstat != '':
                    print(f'Program stopped with message: {runstat}')
                else:
                    print('Program aborted...')
                sys.exit()
            else:
                if zipping:
                    command_zip = (
                        f'gzip -f {star_name}.[l,v,x,y]0000001 '
                        f'{star_name}.b00000 '
                        f'{star_name}.b00001 '
                        f'{star_name}_StrucData_0000001.dat'
                    )
                    os.system(command_zip)
    logfile.close()

    relaunch_advection = [True, 0, 0]

    if calc_dir != '':
        try:
            os.makedirs(calc_dir)
            print(f'Calculation directory {calc_dir} successfully created')
        except OSError:
            print(f'Calculation directory {calc_dir} already exists')
        try:
            os.remove(os.path.join(calc_dir, '*'))
        except OSError:
            pass

    while True:
        with open(
            f'{star_name}.input', 'r', encoding=ENCODING
        ) as input_file:
            inputs = input_file.read()
            ibfile = inputs.rfind(modanfs)+len(modanfs)
            ibfile_end = inputs[ibfile:].find('\n')
            imod = inputs.rfind(nwseqs)+len(nwseqs)
            imod_end = inputs[imod:].find('\n')
            if endphases in inputs:
                requested_stop = True
            if endmodels in inputs:
                requested_stop = True
            iphase = inputs.rfind(phases)+len(phases)
            iphase_end = inputs[iphase:].find('\n')
            modanf = int(inputs[ibfile:ibfile+ibfile_end])
            nwseq = int(inputs[imod:imod+imod_end])
            phase = int(inputs[iphase:iphase+iphase_end])
        if calc_dir != '':
            time_to_transfer = bool(nwseq % ncalc)
        if loop_mode != "" and restart_loop:
            mail_mode = False
            LineLeft = inputs.rfind(deltal)
            LineRight = inputs.rfind(deltat)+14
            NewLine = (
                f'{deltal}{initial_loop[0]:0>7.5f}, '
                f'{deltat}{initial_loop[1]:0>7.5f}'
            )
            to_be_replaced = inputs[LineLeft:LineRight]
            inputs = inputs.replace(to_be_replaced, NewLine)
            with open(f'{star_name}.input', 'w', encoding=ENCODING) as input_file:
                input_file.write(inputs)
            if "up" in loop_mode:
                if initial_loop[1] <= loop_min:
                    initial_loop[0] = initial_loop[0]+loop_step
                    initial_loop[1] = loop_max
                else:
                    initial_loop[1] = initial_loop[1]+loop_step
            if "down" in loop_mode:
                if initial_loop[1] >= loop_max:
                    initial_loop[0] = initial_loop[0]+loop_step
                    initial_loop[1] = loop_min
                else:
                    initial_loop[1] = initial_loop[1]+loop_step

        if phase_stop is not None and phase == phase_stop:
            stop_message = f"{nwseq} : Phase {phase_stop} reached."
            print(stop_message)
            if mail_mode and len(email_receiver) != 0:
                mymail(
                    genec_defaults['email_sender'],
                    email_receiver, stop_message
                )
            elif NOTIFY is True:
                stop_notify(current_dir, stop_message)
            break

        if model_stop is not None and nwseq > model_stop:
            stop_message = f'Model {model_stop} reached.'
            print(stop_message)
            if mail_mode and len(email_receiver) != 0:
                mymail(
                    genec_defaults['email_sender'],
                    email_receiver, stop_message
                )
            elif NOTIFY is True:
                stop_notify(current_dir, stop_message)
            break

        print(f'New model {nwseq} with bfile {modanf}')
        print(f'Phase: {phase}')

        needed_net = [i for i in os.listdir('.') if i[0:3] == 'net']
        needed_for_calc = needed_net+['input_changes.log', f'{star_name}.input']
        if os.path.isfile(f'{star_name}.b{modanf:05d}.gz'):
            needed_for_calc = needed_for_calc+[f'{star_name}.b{modanf:05d}.gz']
        else:
            needed_for_calc = needed_for_calc+[f'{star_name}.b{modanf:05d}']
        try:
            needed_for_calc = needed_for_calc+[f'.PlotData_{star_name}']
        except:
            pass

        if calc_dir != '' and os.getcwd() != calc_dir:
            for file in needed_for_calc:
                shutil.copy2(file, calc_dir)
            os.chdir(calc_dir)

        try:
            os.remove('runfile')
        except OSError:
            pass
        if zipping:
            if os.path.isfile(f'{star_name}.b{modanf:05d}.gz'):
                os.system(f'gunzip {star_name}.b{modanf:05d}.gz')

        os.system(command_launch)
        if zipping:
            command_zip = (
                f'gzip -f '
                f'{star_name}.[l,v,x,y]{nwseq:07d} '
                f'{star_name}.b{modanf:05d} '
                f'{star_name}.b{modanf+1:05d} '
                f'{star_name}_StrucData_{nwseq:07d}.dat'
            )
            os.system(command_zip)

        if calc_dir != '' and time_to_transfer:
            print(f'time to transfer ({nwseq})')
            result_files = [
                i for i in os.listdir('.') if i[0:4] == star_name[0:4]
            ]
            result_files = result_files+['input_changes.log']
            try:
                result_files = result_files+[f'.PlotData_{star_name}']
            except:
                pass
            for file in result_files:
                shutil.move(file, current_dir+'/'+file)
            shutil.copy2('runfile', current_dir)
            os.chdir(current_dir)

        try:
            runlog = open('runfile', 'r', encoding=ENCODING)
        except OSError:
            if (
                "up" in loop_mode
                and (
                    initial_loop[0] >= loop_min
                    or initial_loop[1] >= loop_min
                )
            ):
                restart_loop = True
            elif (
                "down" in loop_mode
                and (
                    initial_loop[0] <= loop_max
                    or initial_loop[1] <= loop_max
                )
            ):
                restart_loop = True
            else:
                stop_message = f'Program stopped with message: {runstat}'
                if mail_mode and len(email_receiver) != 0:
                    mymail(
                        genec_defaults['email_sender'],
                        email_receiver, runstat
                    )
                elif NOTIFY is True:
                    stop_notify(current_dir, stop_message)
                break
            continue
        runstat = runlog.read().strip(' \n\t')
        if runstat != 'running':
            if runstat != '':
                if 'ZAMS reached' in runstat:
                    with open(
                        f'{star_name}.input', 'r', encoding=ENCODING
                    ) as input_file:
                        input_card = input_file.read()
                    nwseq = int(
                        input_card[
                            input_card.rfind(nwseqs)+len(nwseqs):
                            input_card.find('\n modanf')
                        ]
                    )
                    modanf = int(
                        input_card[
                            input_card.rfind(modanfs)+len(modanfs):
                            input_card.find('\n nzmod')
                        ]
                    )
                    try:
                        nzmod = int(
                            input_card[
                                input_card.rfind(nzmods)+len(nzmods):
                                input_card.find('\n&END')
                            ]
                        )
                    except ValueError:
                        nzmod = int(
                            input_card[
                                input_card.rfind(nzmods)+len(nzmods):
                                input_card.find('\n end')
                            ]
                        )
                    print(
                        f'ZAMS reached: NWSEQ= {nwseq}, '
                        f'MODANF= {modanf}, '
                        f'NZMOD= {nzmod}'
                    )
                    if zipping:
                        if os.path.isfile(f'{star_name}.b{modanf:05d}.gz'):
                            os.system(
                                f'gunzip {star_name}.b{modanf:05d}.gz'
                            )
                    os.system(command_launch)
                    if zipping:
                        command_zip = (
                            f'gzip -f '
                            f'{star_name}.[l,v,x,y]{nwseq:07d} '
                            f'{star_name}.b{modanf:05d} '
                            f'{star_name}_StrucData_{nwseq:07d}.dat'
                        )
                        os.system(command_zip)
                    with open('runfile', 'r', encoding=ENCODING) as newrun:
                        runlog = newrun.read().strip(' \n\t')
                        if runlog == 'running':
                            with open(
                                f'{star_name}.input', 'r', encoding=ENCODING
                            ) as input_file:
                                input_card = input_file.read()

                            try:
                                nzmod = int(
                                    input_card[
                                        input_card.rfind(nzmods)
                                        + len(nzmods):input_card.find('\n&END')
                                    ]
                                )
                            except ValueError:
                                nzmod = int(
                                    input_card[
                                        input_card.rfind(nzmods)
                                        + len(nzmods):input_card.find('\n end')
                                    ]
                                )
                            input_card = input_card.replace(
                                f'nzmod={nzmod}', 'nzmod=10')
                            input_card = input_card.replace(
                                'gkorm=.300', 'gkorm=.100')
                            with open(
                                f'{star_name}.input', 'w', encoding=ENCODING
                            ) as input_file:
                                input_file.write(input_card)

                        else:
                            break
                elif 'phase: ' in runstat and requested_stop:
                    stop_message = 'Program reached phase/model requested'
                    if mail_mode and len(email_receiver) != 0:
                        mymail(
                            genec_defaults['email_sender'],
                            email_receiver, runstat)
                    elif NOTIFY is True:
                        stop_notify(current_dir, runstat)
                    break
                elif (
                    'Problem during advection' in runstat
                    or 'Advection not applied' in runstat
                    or 'Problem with conservation of angular momentum during advection' in runstat
                    or 'Ang. mom. variation too large during diffusion' in runstat
                ) and not force_mode:
                    timestep = int(runstat[0:runstat.find(':')])
                    if relaunch_advection[2] == timestep:
                        stop_message = f'Program stopped with message: {runstat}'
                        print('Automatic relaunch for advection failed at previous sequence already twice. ')
                        print(' Please retry with a smaller time step.')
                        if (
                            mail_mode
                            and len(email_receiver) != 0
                        ):
                            mymail(
                                genec_defaults['email_sender'],
                                email_receiver, runstat)
                        elif NOTIFY is True:
                            stop_notify(current_dir, runstat)
                        break
                    relaunch_advection[2] = timestep
                    if relaunch_advection[0]:
                        with open(
                            f'{star_name}.input', 'r', encoding=ENCODING
                        ) as input_card:
                            input_card = input_file.read()
                        nwseq = int(
                            input_card[
                                input_card.rfind(nwseqs)+len(nwseqs):
                                input_card.find('\n modanf')
                            ]
                        )
                        if timestep % 10 == 1:
                            make_command = f"{make_input} {star_name} {nwseq-10}"
                            os.system(make_command)
                        with open(
                            f'{star_name}.input', 'r', encoding=ENCODING
                        ) as input_file:
                            input_card = input_file.read()

                        to_replace = "xcn="+input_card[input_card.rfind(xcns)+len(xcns):input_card.find('\n&END', input_card.rfind(xcns)+len(xcns))]
                        xcn = float(input_card[input_card.rfind(xcns)+len(xcns):input_card.find('\n&END', input_card.rfind(xcns)+len(xcns))])
                        input_card = input_card.replace(to_replace, 'xcn=0.300')
                        with open(
                            f'{star_name}.input', 'w', encoding=ENCODING
                        ) as input_file:
                            input_file.write(input_card)
                        relaunch_advection[1] = relaunch_advection[1] + 1
                        if timestep % 10 == 1 or relaunch_advection[1] > 1:
                            relaunch_advection[0] = False
                    else:
                        stop_message = f'Program stopped with message: {runstat}'
                        if (
                            mail_mode
                            and len(email_receiver) != 0
                        ):
                            mymail(
                                genec_defaults['email_sender'],
                                email_receiver, runstat)
                        elif NOTIFY is True:
                            stop_notify(current_dir, runstat)
                        break
                else:
                    if (
                        "up" in loop_mode
                        and (
                            initial_loop[0] >= loop_min
                            or initial_loop[1] >= loop_min
                        )
                    ):
                        restart_loop = True
                    elif (
                        "down" in loop_mode
                        and (
                            initial_loop[0] <= loop_max
                            or initial_loop[1] <= loop_max
                        )
                    ):
                        restart_loop = True
                    else:
                        stop_message = f'Program stopped with message: {runstat}'
                        if (
                            mail_mode
                            and len(email_receiver) != 0
                        ):
                            mymail(
                                genec_defaults['email_sender'],
                                email_receiver, runstat)
                        elif NOTIFY is True:
                            stop_notify(current_dir, runstat)
                        break
            elif force_mode:
                if (
                    "up" in loop_mode
                    and (
                        initial_loop[0] >= loop_min
                        or initial_loop[1] >= loop_min
                    )
                ):
                    restart_loop = True
                elif (
                    "down" in loop_mode
                    and (
                        initial_loop[0] <= loop_max
                        or initial_loop[1] <= loop_max
                    )
                ):
                    restart_loop = True
                else:
                    stop_message = f'Program stopped with message: {runstat}'
                    if (
                        mail_mode
                        and len(email_receiver) != 0
                    ):
                        mymail(
                            genec_defaults['email_sender'],
                            email_receiver, runstat)
                    elif NOTIFY is True:
                        stop_notify(current_dir, runstat)
                    break
                continue
            else:
                stop_message = 'Program aborted...'
                break

        else:
            restart_loop = False
            mail_mode = base_mail_mode
            if "up" in loop_mode:
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
        result_files = [i for i in os.listdir('.') if i[0:4] == star_name[0:4]]
        result_files = result_files+['input_changes.log']
        try:
            result_files = result_files+[f'.PlotData_{star_name}']
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
    with open('computation.log', 'a', encoding=ENCODING) as logfile:
        logfile.write(
            f'{stop_message}\n'
            f'Computation lasted {computation_time}\n'
            f'----------------------------------------\n\n'
        )

    try:
        os.system('rm buffer_save.dat')
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    main()
