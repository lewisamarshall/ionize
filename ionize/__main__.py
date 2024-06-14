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


@cli.command()
@click.option('--serialize', '-s', is_flag=True)
@click.argument("components", nargs=-1)
def solution(components, serialize):
    ions = components[::2]
    concentrations = [float(c) for c in components[1::2]]
    assert len(ions) == len(concentrations), "There must be a concentration for each ion."
    sol = Solution(ions, concentrations)
    if serialize:
        click.echo(sol.serialize(nested=False, compact=True))
    else:
        click.echo('pH: {:0.4f}'.format(sol.pH))
        click.echo('ionic strength: {:0.4f} M'.format(sol.ionic_strength))
        click.echo('conductivity: {:0.4f} s/M'.format(sol.conductivity()))
        click.echo('Debye length: {:0.4f} um'.format(sol.debye()*1e6))

if __name__ == '__main__':
    cli()
