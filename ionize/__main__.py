import click
from .Solution import Solution
from .Database import Database


@click.group()
def cli():
    pass


@cli.command()
def database():
    click.echo(Database().serialize())


@cli.command()
@click.argument('name')
def ion(name):
    new_ion = Database()[name]
    click.echo(Database()[name].serialize(nested=False, compact=True))


# @cli.command()
# @click.option('--ion', '-i', multiple=True)
# @click.option('--concentration', '-c', multiple=True, type=float)
# def solution(ion, concentration):
#     sol = Solution(ion, concentration)
#     click.echo('pH: {:0.4f}'.format(sol.pH))
#     click.echo('ionic strength: {:0.4f} M'.format(sol.ionic_strength))
#     click.echo('conductivity: {:0.4f} s/M'.format(sol.conductivity()))
#     click.echo('Debye length: {:0.4f} um'.format(sol.debye()*1e6))
#     click.echo(dir(sol))
#
#
# @cli.command()
# def io():
#     pass

if __name__ == '__main__':
    cli()
