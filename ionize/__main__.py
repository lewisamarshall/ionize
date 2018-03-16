import click
from .Solution import Solution
from .Database import Database


@click.group()
def cli():
    pass


@cli.command()
def database():
    '''
    Print the entire ion database.
    '''
    click.echo(Database().serialize())


@cli.command()
@click.argument('name')
def ion(name):
    '''
    Print the physical properties of an ion.
    '''
    new_ion = Database()[name]
    click.echo(Database()[name].serialize(nested = False, compact = True))


@cli.command()
@click.argument('ions')
@click.argument('concentrations')
@click.option('--titrate', '-t', type = (str, float), help = 'titrate the solution with the specified titrant to the specified value of a property, e.g. \'hydrochloric acid\' 8.5')
@click.option('--titration_property', '-p', default = 'pH', help = 'physical property to titrate, default pH')
def solution(ions, concentrations, titrate, titration_property):
    '''
    Compute the physical properties of a solution.
    
    IONS: comma-separated list of ion names, e.g. 'tris,hydrochloric acid'
    CONCENTRATIONS: comma-separated list of concentrations (molar) in same order as ions, e.g. '0.01,0.006'
    '''
    
    initial_solution = Solution(ions.split(','), list(map(float, concentrations.split(','))))
    
    if titrate: click.echo('starting solution:')
    click.echo(initial_solution)
    click.echo(initial_solution.serialize(nested = False, compact = True))
    
    if titrate:
        titrant, target = titrate
        titrated_solution = initial_solution.titrate(titrant, target, titration_property)
        click.echo('\ntitrated solution:')
        click.echo(titrated_solution)
        click.echo(titrated_solution.serialize(nested = False, compact = True))
        click.echo('\ntitration: %f M %s' % ((titrated_solution - initial_solution).concentration(titrant), titrant))


if __name__ == '__main__':
    cli()

