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
    db = Database()
    if name in db:
        new_ion = db[name]
    else:
        search_results = db.search(name)
        if not search_results:
            click.echo('\'%s\' not found' % name)
            return
        click.echo('Did you mean one of these?\n')
        for i in range(len(search_results)):
            click.echo('%i: %s' % (i + 1, search_results[i]))
        which_result = click.prompt('\nEnter the number of the correct ion', type = int)
        new_ion = db[search_results[which_result - 1]]
        click.echo()
    
    click.echo(new_ion.serialize(nested = False, compact = True))


@cli.command()
@click.argument('ions')
@click.argument('concentrations')
@click.option('--titrate', '-t', type = (str, float), default = (None, None), help = 'titrate the solution with the specified titrant to the specified value of a property, e.g. \'hydrochloric acid\' 8.5')
@click.option('--titration_property', '-p', default = 'pH', help = 'physical property to titrate, default pH')
def solution(ions, concentrations, titrate, titration_property):
    '''
    Compute the physical properties of a solution.
    
    IONS: comma-separated list of ion names, e.g. 'tris,hydrochloric acid'
    CONCENTRATIONS: comma-separated list of concentrations (molar) in same order as ions, e.g. '0.01,0.006'
    '''
    
    initial_solution = Solution(ions.split(','), list(map(float, concentrations.split(','))))
    
    if titrate != (None, None): click.echo('starting solution:')
    click.echo(initial_solution)
    click.echo(initial_solution.serialize(nested = False, compact = True))
    
    if titrate != (None, None):
        titrant, target = titrate
        titrated_solution = initial_solution.titrate(titrant, target, titration_property)
        click.echo('\ntitrated solution:')
        click.echo(titrated_solution)
        click.echo(titrated_solution.serialize(nested = False, compact = True))
        click.echo('\ntitration: add %f M %s' % ((titrated_solution - initial_solution).concentration(titrant), titrant))


if __name__ == '__main__':
    cli()

