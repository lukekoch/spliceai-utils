nextflow.enable.dsl=2

params.joblist = "commands.txt"
params.queue = "batch"
params.memory = "4GB"
params.cpus = 1

process runCommand {
    tag "Running command: ${command.take(50)}..."

    input:
    val command

    memory params.memory
    queue params.queue
    cpus params.cpus

    script:
    """
    ${command}
    """
}

workflow {
    def commands = file(params.joblist).readLines()

    channel.from(commands) | runCommand
}